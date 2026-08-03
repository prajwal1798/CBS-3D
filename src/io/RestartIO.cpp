//=============================================================================
// CBS3D++_SI
//
// Distributed native restart/checkpoint implementation and one-time importer
// for the existing ASCII VTU pieces.
//=============================================================================

#include "cbs/io/RestartIO.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    namespace
    {
        namespace fs = std::filesystem;

        constexpr std::uint32_t restart_version = 1U;
        constexpr std::uint64_t restart_end_marker = 0xC0DEC0DEC0DEC0DEULL;

        constexpr std::array<char, 16> restart_magic =
        {
            'C', 'B', 'S', '3', 'D', '_', 'R', 'E',
            'S', 'T', 'A', 'R', 'T', '_', '1', '\0'
        };

        bool environment_flag_enabled(const char* name)
        {
            const char* value = std::getenv(name);

            return value != nullptr &&
                   value[0] != '\0' &&
                   std::string(value) != "0";
        }

        std::string environment_value(
            const char* name,
            const std::string& fallback = "")
        {
            const char* value = std::getenv(name);

            if (value == nullptr || value[0] == '\0')
            {
                return fallback;
            }

            return std::string(value);
        }

        Int parse_nonnegative_environment_int(
            const char* name,
            const Int fallback)
        {
            const std::string text = environment_value(name);

            if (text.empty())
            {
                return fallback;
            }

            std::size_t parsed = 0;
            long long value = 0;

            try
            {
                value = std::stoll(text, &parsed, 10);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error(
                    std::string("Invalid integer in environment variable ") + name);
            }

            if (parsed != text.size() ||
                value < 0 ||
                value > std::numeric_limits<Int>::max())
            {
                throw std::runtime_error(
                    std::string("Environment variable out of range: ") + name);
            }

            return static_cast<Int>(value);
        }

        std::string lower_copy(std::string value)
        {
            for (char& character : value)
            {
                character = static_cast<char>(
                    std::tolower(static_cast<unsigned char>(character)));
            }

            return value;
        }

        std::string case_name_from_local_path(
            const std::string& rank_local_case_name)
        {
            const fs::path local_path(rank_local_case_name);
            std::string name = local_path.filename().string();

            const std::string marker = "_rank_";
            const std::size_t marker_position = name.rfind(marker);

            if (marker_position != std::string::npos)
            {
                const std::string suffix =
                    name.substr(marker_position + marker.size());

                const bool numeric_suffix =
                    !suffix.empty() &&
                    std::all_of(
                        suffix.begin(),
                        suffix.end(),
                        [](const unsigned char character)
                        {
                            return character >= '0' && character <= '9';
                        });

                if (numeric_suffix)
                {
                    name.resize(marker_position);
                }
            }

            if (name.empty())
            {
                throw std::runtime_error(
                    "RestartIO obtained an empty distributed case name");
            }

            return name;
        }

        std::string step_tag(const Int iteration)
        {
            std::ostringstream text;
            text << "step_" << std::setw(8) << std::setfill('0') << iteration;
            return text.str();
        }

        std::string rank_tag(const Int rank)
        {
            std::ostringstream text;
            text << "rank_" << std::setw(4) << std::setfill('0') << rank;
            return text.str();
        }

        std::string native_rank_file_name(
            const std::string& case_name,
            const Int iteration,
            const Int rank)
        {
            return case_name + "_" + step_tag(iteration)
                 + "_" + rank_tag(rank) + ".cbsrst";
        }

        std::string legacy_vtu_file_name(
            const std::string& case_name,
            const Int iteration,
            const Int rank)
        {
            return case_name + "_" + step_tag(iteration)
                 + "_" + rank_tag(rank) + ".vtu";
        }

        std::uint64_t local_node_map_hash(const CBSStateSI& s)
        {
            // FNV-1a over the exact rank-local node ordering.  A restart is
            // rejected when the partition exporter or local numbering changed.
            std::uint64_t hash = 1469598103934665603ULL;

            const auto mix_byte = [&hash](const std::uint8_t byte)
            {
                hash ^= byte;
                hash *= 1099511628211ULL;
            };

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const std::size_t index = static_cast<std::size_t>(ip);

                if (index >= s.local_to_global_node.size())
                {
                    throw std::runtime_error(
                        "RestartIO found an incomplete local-to-global node map");
                }

                const std::uint64_t value = static_cast<std::uint64_t>(
                    s.local_to_global_node[index]);

                for (int byte = 0; byte < 8; ++byte)
                {
                    mix_byte(static_cast<std::uint8_t>(
                        (value >> (8 * byte)) & 0xFFU));
                }
            }

            return hash;
        }

        template <typename T>
        void write_binary_value(std::ostream& output, const T& value)
        {
            output.write(
                reinterpret_cast<const char*>(&value),
                static_cast<std::streamsize>(sizeof(T)));

            if (!output)
            {
                throw std::runtime_error(
                    "RestartIO failed while writing a binary value");
            }
        }

        template <typename T>
        T read_binary_value(std::istream& input)
        {
            T value{};

            input.read(
                reinterpret_cast<char*>(&value),
                static_cast<std::streamsize>(sizeof(T)));

            if (!input)
            {
                throw std::runtime_error(
                    "RestartIO encountered a truncated binary restart file");
            }

            return value;
        }

        void write_native_rank_file(
            const CBSStateSI& s,
            const fs::path& path,
            const Int completed_iteration)
        {
            std::ofstream output(
                path,
                std::ios::binary | std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "RestartIO cannot create rank-local checkpoint: "
                    + path.string());
            }

            output.write(
                restart_magic.data(),
                static_cast<std::streamsize>(restart_magic.size()));

            if (!output)
            {
                throw std::runtime_error(
                    "RestartIO failed while writing checkpoint magic");
            }

            const std::int32_t rank = s.mpi_rank;
            const std::int32_t mpi_size = s.mpi_size;
            const std::int64_t local_npoin = s.cfg.npoin;
            const std::int64_t global_npoin =
                s.partition_metadata.global_npoin;
            const std::int64_t iteration = completed_iteration;
            const std::uint64_t map_hash = local_node_map_hash(s);

            std::uint32_t flags = 0U;

            if (s.cfg.temp_calc > 0)
            {
                flags |= 1U;
            }

            if (s.cfg.turbulence_on > 0)
            {
                flags |= 2U;
            }

            write_binary_value(output, restart_version);
            write_binary_value(output, rank);
            write_binary_value(output, mpi_size);
            write_binary_value(output, local_npoin);
            write_binary_value(output, global_npoin);
            write_binary_value(output, iteration);
            write_binary_value(output, s.cfg.rtime);
            write_binary_value(output, s.cfg.dtreal);
            write_binary_value(output, map_hash);
            write_binary_value(output, flags);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const std::int64_t global_id =
                    s.local_to_global_node[static_cast<std::size_t>(ip)];

                write_binary_value(output, global_id);
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                for (Int dimension = 1;
                     dimension <= s.cfg.ndim;
                     ++dimension)
                {
                    write_binary_value(output, s.unkno(dimension, ip));
                }
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                write_binary_value(output, s.pres(ip));
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                write_binary_value(output, s.temperature(ip));
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                write_binary_value(output, s.nu_tilde(ip));
            }

            write_binary_value(output, restart_end_marker);

            output.flush();

            if (!output)
            {
                throw std::runtime_error(
                    "RestartIO failed while flushing rank-local checkpoint");
            }
        }

        struct Manifest
        {
            std::string case_name;
            Int completed_iteration = 0;
            Int mpi_size = 0;
            std::int64_t global_npoin = 0;
            std::int64_t global_nelem = 0;
            Real physical_time = 0.0;
            Real time_step = 0.0;
        };

        Manifest read_manifest(const fs::path& checkpoint_directory)
        {
            const fs::path manifest_path =
                checkpoint_directory / "restart_manifest.txt";

            std::ifstream input(manifest_path);

            if (!input)
            {
                throw std::runtime_error(
                    "RestartIO cannot open restart manifest: "
                    + manifest_path.string());
            }

            Manifest manifest;
            std::string key;

            while (input >> key)
            {
                if (key == "format_version")
                {
                    std::uint32_t version = 0;
                    input >> version;

                    if (version != restart_version)
                    {
                        throw std::runtime_error(
                            "RestartIO manifest format version is unsupported");
                    }
                }
                else if (key == "case_name")
                {
                    input >> manifest.case_name;
                }
                else if (key == "completed_iteration")
                {
                    input >> manifest.completed_iteration;
                }
                else if (key == "mpi_size")
                {
                    input >> manifest.mpi_size;
                }
                else if (key == "global_npoin")
                {
                    input >> manifest.global_npoin;
                }
                else if (key == "global_nelem")
                {
                    input >> manifest.global_nelem;
                }
                else if (key == "physical_time")
                {
                    input >> manifest.physical_time;
                }
                else if (key == "time_step")
                {
                    input >> manifest.time_step;
                }
                else
                {
                    std::string ignored;
                    std::getline(input, ignored);
                }

                if (!input)
                {
                    throw std::runtime_error(
                        "RestartIO found malformed restart manifest data");
                }
            }

            if (manifest.case_name.empty() ||
                manifest.completed_iteration < 0 ||
                manifest.mpi_size < 1 ||
                manifest.global_npoin < 1 ||
                manifest.global_nelem < 1)
            {
                throw std::runtime_error(
                    "RestartIO restart manifest is incomplete");
            }

            return manifest;
        }

        RestartIO::LoadResult read_native_checkpoint(
            CBSStateSI& s,
            const std::string& rank_local_case_name,
            const fs::path& checkpoint_directory)
        {
            const std::string case_name =
                case_name_from_local_path(rank_local_case_name);

            const Manifest manifest = read_manifest(checkpoint_directory);

            if (manifest.case_name != case_name)
            {
                throw std::runtime_error(
                    "RestartIO case name does not match restart manifest");
            }

            if (manifest.mpi_size != s.mpi_size)
            {
                throw std::runtime_error(
                    "RestartIO requires the same MPI process count as the checkpoint");
            }

            if (manifest.global_npoin != s.partition_metadata.global_npoin ||
                manifest.global_nelem != s.partition_metadata.global_nelem)
            {
                throw std::runtime_error(
                    "RestartIO global mesh size does not match the checkpoint");
            }

            const fs::path rank_path =
                checkpoint_directory /
                native_rank_file_name(
                    case_name,
                    manifest.completed_iteration,
                    s.mpi_rank);

            std::ifstream input(rank_path, std::ios::binary);

            if (!input)
            {
                throw std::runtime_error(
                    "RestartIO cannot open rank-local checkpoint: "
                    + rank_path.string());
            }

            std::array<char, 16> magic{};
            input.read(
                magic.data(),
                static_cast<std::streamsize>(magic.size()));

            if (!input || magic != restart_magic)
            {
                throw std::runtime_error(
                    "RestartIO rank file has invalid magic");
            }

            const std::uint32_t version =
                read_binary_value<std::uint32_t>(input);
            const std::int32_t rank =
                read_binary_value<std::int32_t>(input);
            const std::int32_t mpi_size =
                read_binary_value<std::int32_t>(input);
            const std::int64_t local_npoin =
                read_binary_value<std::int64_t>(input);
            const std::int64_t global_npoin =
                read_binary_value<std::int64_t>(input);
            const std::int64_t completed_iteration =
                read_binary_value<std::int64_t>(input);
            const Real physical_time =
                read_binary_value<Real>(input);
            const Real time_step =
                read_binary_value<Real>(input);
            const std::uint64_t stored_map_hash =
                read_binary_value<std::uint64_t>(input);
            const std::uint32_t flags =
                read_binary_value<std::uint32_t>(input);

            if (version != restart_version ||
                rank != s.mpi_rank ||
                mpi_size != s.mpi_size ||
                local_npoin != s.cfg.npoin ||
                global_npoin != s.partition_metadata.global_npoin ||
                completed_iteration != manifest.completed_iteration)
            {
                throw std::runtime_error(
                    "RestartIO rank-file metadata does not match the current partition");
            }

            if (stored_map_hash != local_node_map_hash(s))
            {
                throw std::runtime_error(
                    "RestartIO local node ordering differs from the checkpoint");
            }

            const bool checkpoint_temperature_enabled = (flags & 1U) != 0U;
            const bool checkpoint_turbulence_enabled = (flags & 2U) != 0U;

            if (checkpoint_temperature_enabled != (s.cfg.temp_calc > 0))
            {
                throw std::runtime_error(
                    "RestartIO temperature-equation setting differs from the checkpoint");
            }

            if (checkpoint_turbulence_enabled != (s.cfg.turbulence_on > 0))
            {
                throw std::runtime_error(
                    "RestartIO turbulence setting differs from the checkpoint");
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const std::int64_t stored_global_id =
                    read_binary_value<std::int64_t>(input);

                if (stored_global_id !=
                    s.local_to_global_node[static_cast<std::size_t>(ip)])
                {
                    throw std::runtime_error(
                        "RestartIO local/global node map differs from checkpoint");
                }
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                for (Int dimension = 1;
                     dimension <= s.cfg.ndim;
                     ++dimension)
                {
                    s.unkno(dimension, ip) =
                        read_binary_value<Real>(input);
                }
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.pres(ip) = read_binary_value<Real>(input);
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.temperature(ip) = read_binary_value<Real>(input);
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.nu_tilde(ip) = read_binary_value<Real>(input);
            }

            const std::uint64_t end_marker =
                read_binary_value<std::uint64_t>(input);

            if (end_marker != restart_end_marker)
            {
                throw std::runtime_error(
                    "RestartIO rank file is truncated or corrupt");
            }

            RestartIO::LoadResult result;
            result.loaded = true;
            result.imported_from_legacy_vtu = false;
            result.completed_iteration =
                static_cast<Int>(completed_iteration);
            result.physical_time = physical_time;
            result.time_step = time_step;
            return result;
        }

        std::string read_text_file(const fs::path& path)
        {
            std::ifstream input(path, std::ios::binary);

            if (!input)
            {
                throw std::runtime_error(
                    "RestartIO cannot open legacy VTU piece: " + path.string());
            }

            std::ostringstream buffer;
            buffer << input.rdbuf();

            if (!input.good() && !input.eof())
            {
                throw std::runtime_error(
                    "RestartIO failed while reading legacy VTU piece");
            }

            return buffer.str();
        }

        std::string data_array_payload(
            const std::string& xml,
            const std::string& name)
        {
            const std::string marker = "Name=\"" + name + "\"";
            const std::size_t marker_position = xml.find(marker);

            if (marker_position == std::string::npos)
            {
                throw std::runtime_error(
                    "RestartIO legacy VTU is missing DataArray " + name);
            }

            const std::size_t payload_begin = xml.find('>', marker_position);

            if (payload_begin == std::string::npos)
            {
                throw std::runtime_error(
                    "RestartIO found malformed legacy VTU DataArray " + name);
            }

            const std::size_t payload_end =
                xml.find("</DataArray>", payload_begin + 1U);

            if (payload_end == std::string::npos)
            {
                throw std::runtime_error(
                    "RestartIO found unterminated legacy VTU DataArray " + name);
            }

            return xml.substr(
                payload_begin + 1U,
                payload_end - payload_begin - 1U);
        }

        std::vector<Real> parse_real_array(
            const std::string& payload,
            const std::size_t expected_values,
            const std::string& field_name)
        {
            std::istringstream input(payload);
            std::vector<Real> values(expected_values, 0.0);

            for (std::size_t index = 0;
                 index < expected_values;
                 ++index)
            {
                if (!(input >> values[index]) ||
                    !std::isfinite(values[index]))
                {
                    throw std::runtime_error(
                        "RestartIO legacy VTU contains invalid " + field_name);
                }
            }

            Real extra = 0.0;

            if (input >> extra)
            {
                throw std::runtime_error(
                    "RestartIO legacy VTU contains too many values in "
                    + field_name);
            }

            return values;
        }

        std::vector<std::int64_t> parse_integer_array(
            const std::string& payload,
            const std::size_t expected_values,
            const std::string& field_name)
        {
            std::istringstream input(payload);
            std::vector<std::int64_t> values(expected_values, 0);

            for (std::size_t index = 0;
                 index < expected_values;
                 ++index)
            {
                if (!(input >> values[index]))
                {
                    throw std::runtime_error(
                        "RestartIO legacy VTU contains invalid " + field_name);
                }
            }

            std::int64_t extra = 0;

            if (input >> extra)
            {
                throw std::runtime_error(
                    "RestartIO legacy VTU contains too many values in "
                    + field_name);
            }

            return values;
        }

        RestartIO::LoadResult read_legacy_vtu_checkpoint(
            CBSStateSI& s,
            const std::string& rank_local_case_name,
            const fs::path& legacy_directory,
            const Int completed_iteration)
        {
            if (s.cfg.transient_on > 0)
            {
                throw std::runtime_error(
                    "RestartIO legacy VTU import is allowed only for the current "
                    "steady pseudo-time calculation; transient history is absent");
            }

            const std::string case_name =
                case_name_from_local_path(rank_local_case_name);

            const fs::path piece_path =
                legacy_directory /
                legacy_vtu_file_name(
                    case_name,
                    completed_iteration,
                    s.mpi_rank);

            const std::string xml = read_text_file(piece_path);
            const std::size_t npoin = static_cast<std::size_t>(s.cfg.npoin);

            const std::vector<std::int64_t> global_ids =
                parse_integer_array(
                    data_array_payload(xml, "global_node_id"),
                    npoin,
                    "global_node_id");

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (global_ids[static_cast<std::size_t>(ip - 1)] !=
                    s.local_to_global_node[static_cast<std::size_t>(ip)])
                {
                    throw std::runtime_error(
                        "RestartIO legacy VTU node ordering differs from the "
                        "current rank-local partition");
                }
            }

            const std::vector<Real> velocity =
                parse_real_array(
                    data_array_payload(xml, "velocity"),
                    3U * npoin,
                    "velocity");

            const std::vector<Real> pressure =
                parse_real_array(
                    data_array_payload(xml, "pressure"),
                    npoin,
                    "pressure");

            const std::vector<Real> temperature =
                parse_real_array(
                    data_array_payload(xml, "temperature"),
                    npoin,
                    "temperature");

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const std::size_t node = static_cast<std::size_t>(ip - 1);

                s.unkno(1, ip) = velocity[3U * node];
                s.unkno(2, ip) = velocity[3U * node + 1U];
                s.unkno(3, ip) = velocity[3U * node + 2U];
                s.pres(ip) = pressure[node];
                s.temperature(ip) = temperature[node];
            }

            s.nu_tilde.fill(0.0);

            RestartIO::LoadResult result;
            result.loaded = true;
            result.imported_from_legacy_vtu = true;
            result.completed_iteration = completed_iteration;
            result.physical_time = 0.0;
            result.time_step = s.cfg.dtreal;
            return result;
        }

        void initialise_restart_histories(
            CBSStateSI& s,
            const RestartIO::LoadResult& result)
        {
            s.unkn1 = s.unkno;
            s.unknn1 = s.unkno;
            s.unknn2 = s.unkno;

            s.pres1 = s.pres;

            s.temperature1 = s.temperature;
            s.tempert1 = s.temperature;
            s.tempert2 = s.temperature;

            s.nu_tilde1 = s.nu_tilde;

            s.cfg.rtime = result.physical_time;
            s.cfg.dtreal = result.time_step;
            s.cfg.iiter_total = result.completed_iteration;
            s.cfg.iitime_start = result.completed_iteration + 1;
            s.cfg.istart = 2;
        }

#ifdef CBS3D_USE_MPI
        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("RestartIO MPI failure in ") + operation);
            }
        }
#endif
    }

    RestartIO::LoadResult RestartIO::loadIfRequested(
        CBSStateSI& s,
        const std::string& rank_local_case_name)
    {
        const bool requested =
            s.cfg.restart_opt > 0 ||
            environment_flag_enabled("CBS3D_RESTART");

        if (!requested)
        {
            return LoadResult{};
        }

        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "RestartIO distributed restart requires more than one MPI rank");
        }

        const std::string root_text =
            environment_value("CBS3D_RESTART_ROOT");

        if (root_text.empty())
        {
            throw std::runtime_error(
                "RestartIO requires CBS3D_RESTART_ROOT when restart is enabled");
        }

        const std::string format = lower_copy(
            environment_value("CBS3D_RESTART_FORMAT", "native"));

        LoadResult result;

        if (format == "native")
        {
            result = read_native_checkpoint(
                s,
                rank_local_case_name,
                fs::path(root_text));
        }
        else if (format == "legacy_vtu" ||
                 format == "vtu" ||
                 format == "legacy")
        {
            const Int completed_iteration =
                parse_nonnegative_environment_int(
                    "CBS3D_RESTART_ITERATION",
                    -1);

            if (completed_iteration < 1)
            {
                throw std::runtime_error(
                    "RestartIO legacy VTU import requires "
                    "CBS3D_RESTART_ITERATION >= 1");
            }

            result = read_legacy_vtu_checkpoint(
                s,
                rank_local_case_name,
                fs::path(root_text),
                completed_iteration);
        }
        else
        {
            throw std::runtime_error(
                "RestartIO CBS3D_RESTART_FORMAT must be native or legacy_vtu");
        }

        initialise_restart_histories(s, result);
        return result;
    }

    Int RestartIO::checkpointInterval()
    {
        return parse_nonnegative_environment_int(
            "CBS3D_CHECKPOINT_EVERY",
            0);
    }

    std::string RestartIO::checkpointRoot(
        const std::string& rank_local_case_name)
    {
        const std::string configured =
            environment_value("CBS3D_CHECKPOINT_ROOT");

        if (!configured.empty())
        {
            return configured;
        }

        const std::string case_name =
            case_name_from_local_path(rank_local_case_name);

        return (fs::path("output") / case_name / "restart").string();
    }

    void RestartIO::writeCheckpoint(
        const CBSStateSI& s,
        const std::string& rank_local_case_name,
        const Int completed_iteration)
    {
#ifdef CBS3D_USE_MPI
        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "RestartIO::writeCheckpoint requires more than one MPI rank");
        }

        if (completed_iteration < 0)
        {
            throw std::runtime_error(
                "RestartIO cannot write a negative completed iteration");
        }

        const std::string case_name =
            case_name_from_local_path(rank_local_case_name);

        const fs::path checkpoint_directory =
            fs::path(checkpointRoot(rank_local_case_name)) /
            step_tag(completed_iteration);

        int directory_ok = 1;

        if (s.mpi_rank == 0)
        {
            try
            {
                fs::create_directories(checkpoint_directory);
            }
            catch (const std::exception&)
            {
                directory_ok = 0;
            }
        }

        check_mpi(
            MPI_Bcast(
                &directory_ok,
                1,
                MPI_INT,
                0,
                MPI_COMM_WORLD),
            "MPI_Bcast checkpoint-directory status");

        if (directory_ok == 0)
        {
            throw std::runtime_error(
                "RestartIO failed to create checkpoint directory");
        }

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier checkpoint directory");

        const fs::path final_rank_path =
            checkpoint_directory /
            native_rank_file_name(
                case_name,
                completed_iteration,
                s.mpi_rank);

        const fs::path temporary_rank_path =
            final_rank_path.string() + ".tmp";

        int local_write_ok = 1;
        std::string local_write_error;

        try
        {
            write_native_rank_file(
                s,
                temporary_rank_path,
                completed_iteration);
        }
        catch (const std::exception& error)
        {
            local_write_ok = 0;
            local_write_error = error.what();
        }

        int global_write_ok = 0;

        check_mpi(
            MPI_Allreduce(
                &local_write_ok,
                &global_write_ok,
                1,
                MPI_INT,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce checkpoint write status");

        if (global_write_ok == 0)
        {
            std::error_code ignored;
            fs::remove(temporary_rank_path, ignored);

            if (!local_write_ok && !local_write_error.empty())
            {
                throw std::runtime_error(local_write_error);
            }

            throw std::runtime_error(
                "RestartIO one or more ranks failed to write checkpoint data");
        }

        int local_rename_ok = 1;

        try
        {
            std::error_code ignored;
            fs::remove(final_rank_path, ignored);
            fs::rename(temporary_rank_path, final_rank_path);
        }
        catch (const std::exception&)
        {
            local_rename_ok = 0;
        }

        int global_rename_ok = 0;

        check_mpi(
            MPI_Allreduce(
                &local_rename_ok,
                &global_rename_ok,
                1,
                MPI_INT,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce checkpoint rename status");

        if (global_rename_ok == 0)
        {
            throw std::runtime_error(
                "RestartIO failed while committing rank checkpoint files");
        }

        int manifest_ok = 1;

        if (s.mpi_rank == 0)
        {
            try
            {
                const fs::path manifest_path =
                    checkpoint_directory / "restart_manifest.txt";
                const fs::path temporary_manifest_path =
                    manifest_path.string() + ".tmp";

                std::ofstream manifest(
                    temporary_manifest_path,
                    std::ios::trunc);

                if (!manifest)
                {
                    throw std::runtime_error(
                        "cannot create checkpoint manifest");
                }

                manifest << std::setprecision(17);
                manifest << "format_version " << restart_version << '\n';
                manifest << "case_name " << case_name << '\n';
                manifest << "completed_iteration "
                         << completed_iteration << '\n';
                manifest << "mpi_size " << s.mpi_size << '\n';
                manifest << "global_npoin "
                         << s.partition_metadata.global_npoin << '\n';
                manifest << "global_nelem "
                         << s.partition_metadata.global_nelem << '\n';
                manifest << "physical_time " << s.cfg.rtime << '\n';
                manifest << "time_step " << s.cfg.dtreal << '\n';
                manifest.flush();

                if (!manifest)
                {
                    throw std::runtime_error(
                        "cannot flush checkpoint manifest");
                }

                manifest.close();

                std::error_code ignored;
                fs::remove(manifest_path, ignored);
                fs::rename(temporary_manifest_path, manifest_path);
            }
            catch (const std::exception&)
            {
                manifest_ok = 0;
            }
        }

        check_mpi(
            MPI_Bcast(
                &manifest_ok,
                1,
                MPI_INT,
                0,
                MPI_COMM_WORLD),
            "MPI_Bcast checkpoint-manifest status");

        if (manifest_ok == 0)
        {
            throw std::runtime_error(
                "RestartIO failed while committing checkpoint manifest");
        }

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier completed checkpoint");
#else
        (void)s;
        (void)rank_local_case_name;
        (void)completed_iteration;
        throw std::runtime_error(
            "RestartIO::writeCheckpoint requires an MPI-enabled build");
#endif
    }
}
