#!/usr/bin/env python3
"""Apply DD-5B1 appended-raw-binary VTU output to DistributedPost.cpp.

The transformation is intentionally narrow:

* numerical assembly and the distributed solver loop are untouched;
* PVTU/PVD filenames and field names remain unchanged;
* only rank-local VTU payload encoding changes from ASCII to VTK XML appended raw binary;
* steady pseudo-time PVD entries use the iteration number instead of repeated zero.

Compatible with Python 3.6 on Swansea Sunbird.
"""

from pathlib import Path
import sys

SOURCE = Path("src/io/DistributedPost.cpp")

INCLUDE_OLD = """#include <algorithm>\n#include <cmath>\n#include <filesystem>\n#include <fstream>\n#include <iomanip>\n#include <sstream>\n#include <stdexcept>\n#include <string>\n"""

INCLUDE_NEW = """#include <algorithm>\n#include <array>\n#include <bit>\n#include <cmath>\n#include <cstdint>\n#include <cstring>\n#include <filesystem>\n#include <fstream>\n#include <iomanip>\n#include <limits>\n#include <sstream>\n#include <stdexcept>\n#include <string>\n#include <type_traits>\n"""

START_MARKER = "        void write_piece(\n"
END_MARKER = "        void write_parallel_descriptor(\n"

REPLACEMENT = r'''        using VtuHeader = std::uint64_t;

        constexpr std::size_t binary_output_buffer_bytes = 1024U * 1024U;

        enum VtuBlock : std::size_t
        {
            points_block = 0,
            connectivity_block,
            offsets_block,
            types_block,
            pressure_block,
            temperature_block,
            velocity_block,
            velocity_magnitude_block,
            global_node_id_block,
            owner_rank_block,
            is_owned_block,
            ghost_type_block,
            node_domain_kind_block,
            velocity_bc_type_block,
            pressure_fixed_block,
            cell_kind_block,
            global_element_id_block,
            material_id_block,
            bc_id_block,
            parent_global_element_block,
            vtu_block_count
        };

        std::uint64_t checked_binary_product(
            const std::uint64_t count,
            const std::uint64_t item_bytes,
            const char* description)
        {
            if (item_bytes != 0 &&
                count > std::numeric_limits<std::uint64_t>::max() / item_bytes)
            {
                throw std::runtime_error(
                    std::string("DistributedPost binary-size overflow for ") +
                    description);
            }

            return count * item_bytes;
        }

        std::uint64_t checked_binary_sum(
            const std::uint64_t left,
            const std::uint64_t right,
            const char* description)
        {
            if (left > std::numeric_limits<std::uint64_t>::max() - right)
            {
                throw std::runtime_error(
                    std::string("DistributedPost binary-offset overflow for ") +
                    description);
            }

            return left + right;
        }

        std::int32_t vtk_local_index(const Int one_based_index)
        {
            if (one_based_index < 1 ||
                static_cast<std::uint64_t>(one_based_index - 1) >
                    static_cast<std::uint64_t>(
                        std::numeric_limits<std::int32_t>::max()))
            {
                throw std::runtime_error(
                    "DistributedPost local VTK index is outside Int32 range");
            }

            return static_cast<std::int32_t>(one_based_index - 1);
        }

        std::int32_t vtk_int32(const Int value, const char* description)
        {
            if (value < std::numeric_limits<std::int32_t>::min() ||
                value > std::numeric_limits<std::int32_t>::max())
            {
                throw std::runtime_error(
                    std::string("DistributedPost ") + description +
                    " is outside Int32 range");
            }

            return static_cast<std::int32_t>(value);
        }

        class AppendedBlockWriter
        {
        public:
            AppendedBlockWriter(
                std::ostream& output,
                const std::uint64_t expected_bytes)
                : output_(output),
                  expected_bytes_(expected_bytes)
            {
                output_.write(
                    reinterpret_cast<const char*>(&expected_bytes_),
                    static_cast<std::streamsize>(sizeof(expected_bytes_)));

                if (!output_)
                {
                    throw std::runtime_error(
                        "DistributedPost cannot write a VTU binary block header");
                }
            }

            template <typename Value>
            void append(const Value& value)
            {
                static_assert(
                    std::is_trivially_copyable_v<Value>,
                    "VTU binary values must be trivially copyable");

                if (buffer_used_ + sizeof(Value) > buffer_.size())
                {
                    flush();
                }

                std::memcpy(
                    buffer_.data() + buffer_used_,
                    &value,
                    sizeof(Value));

                buffer_used_ += sizeof(Value);
                written_bytes_ = checked_binary_sum(
                    written_bytes_,
                    sizeof(Value),
                    "VTU binary block payload");

                if (written_bytes_ > expected_bytes_)
                {
                    throw std::runtime_error(
                        "DistributedPost wrote beyond a VTU binary block size");
                }
            }

            void finish()
            {
                flush();

                if (written_bytes_ != expected_bytes_)
                {
                    throw std::runtime_error(
                        "DistributedPost VTU binary block size mismatch");
                }
            }

        private:
            void flush()
            {
                if (buffer_used_ == 0)
                {
                    return;
                }

                output_.write(
                    buffer_.data(),
                    static_cast<std::streamsize>(buffer_used_));

                if (!output_)
                {
                    throw std::runtime_error(
                        "DistributedPost cannot write a VTU binary block payload");
                }

                buffer_used_ = 0;
            }

            std::ostream& output_;
            std::uint64_t expected_bytes_ = 0;
            std::uint64_t written_bytes_ = 0;
            std::array<char, binary_output_buffer_bytes> buffer_{};
            std::size_t buffer_used_ = 0;
        };

        void write_piece(
            const CBSStateSI& s,
            const std::string& case_name,
            const Int iteration)
        {
            if constexpr (std::endian::native != std::endian::little)
            {
                throw std::runtime_error(
                    "DistributedPost appended raw VTU requires little-endian hardware");
            }

            if (s.cfg.npoin < 0 || s.cfg.nelem < 0 || s.cfg.nboun < 0)
            {
                throw std::runtime_error(
                    "DistributedPost obtained negative local mesh sizes");
            }

            const fs::path output_path =
                case_output_directory(case_name) /
                piece_file_name(case_name, iteration, s.mpi_rank);

            std::ofstream output(
                output_path,
                std::ios::binary | std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost cannot create rank-local VTU piece");
            }

            const std::uint64_t point_count =
                static_cast<std::uint64_t>(s.cfg.npoin);
            const std::uint64_t element_count =
                static_cast<std::uint64_t>(s.cfg.nelem);
            const std::uint64_t boundary_count =
                static_cast<std::uint64_t>(s.cfg.nboun);
            const std::uint64_t total_cells = checked_binary_sum(
                element_count,
                boundary_count,
                "VTU cell count");
            const std::uint64_t connectivity_values = checked_binary_sum(
                checked_binary_product(
                    element_count,
                    4,
                    "tetrahedral connectivity count"),
                checked_binary_product(
                    boundary_count,
                    3,
                    "boundary connectivity count"),
                "VTU connectivity count");

            std::array<std::uint64_t, vtu_block_count> block_sizes{};
            block_sizes[points_block] = checked_binary_product(
                checked_binary_product(point_count, 3, "point components"),
                sizeof(double),
                "point coordinates");
            block_sizes[connectivity_block] = checked_binary_product(
                connectivity_values,
                sizeof(std::int32_t),
                "connectivity");
            block_sizes[offsets_block] = checked_binary_product(
                total_cells,
                sizeof(std::int64_t),
                "cell offsets");
            block_sizes[types_block] = checked_binary_product(
                total_cells,
                sizeof(std::uint8_t),
                "cell types");
            block_sizes[pressure_block] = checked_binary_product(
                point_count,
                sizeof(double),
                "pressure");
            block_sizes[temperature_block] = checked_binary_product(
                point_count,
                sizeof(double),
                "temperature");
            block_sizes[velocity_block] = checked_binary_product(
                checked_binary_product(point_count, 3, "velocity components"),
                sizeof(double),
                "velocity");
            block_sizes[velocity_magnitude_block] = checked_binary_product(
                point_count,
                sizeof(double),
                "velocity magnitude");
            block_sizes[global_node_id_block] = checked_binary_product(
                point_count,
                sizeof(std::int64_t),
                "global node IDs");
            block_sizes[owner_rank_block] = checked_binary_product(
                point_count,
                sizeof(std::int32_t),
                "owner ranks");
            block_sizes[is_owned_block] = checked_binary_product(
                point_count,
                sizeof(std::uint8_t),
                "owned-node flags");
            block_sizes[ghost_type_block] = checked_binary_product(
                point_count,
                sizeof(std::uint8_t),
                "VTK ghost flags");
            block_sizes[node_domain_kind_block] = checked_binary_product(
                point_count,
                sizeof(std::int32_t),
                "node domain kinds");
            block_sizes[velocity_bc_type_block] = checked_binary_product(
                point_count,
                sizeof(std::int32_t),
                "velocity boundary types");
            block_sizes[pressure_fixed_block] = checked_binary_product(
                point_count,
                sizeof(std::uint8_t),
                "fixed-pressure flags");
            block_sizes[cell_kind_block] = checked_binary_product(
                total_cells,
                sizeof(std::uint8_t),
                "cell kinds");
            block_sizes[global_element_id_block] = checked_binary_product(
                total_cells,
                sizeof(std::int64_t),
                "global element IDs");
            block_sizes[material_id_block] = checked_binary_product(
                total_cells,
                sizeof(std::int32_t),
                "material IDs");
            block_sizes[bc_id_block] = checked_binary_product(
                total_cells,
                sizeof(std::int32_t),
                "boundary IDs");
            block_sizes[parent_global_element_block] = checked_binary_product(
                total_cells,
                sizeof(std::int64_t),
                "parent global elements");

            std::array<std::uint64_t, vtu_block_count> block_offsets{};
            std::uint64_t next_offset = 0;

            for (std::size_t block = 0; block < block_sizes.size(); ++block)
            {
                block_offsets[block] = next_offset;
                next_offset = checked_binary_sum(
                    next_offset,
                    checked_binary_sum(
                        sizeof(VtuHeader),
                        block_sizes[block],
                        "VTU block extent"),
                    "VTU appended offset");
            }

            output << "<?xml version=\"1.0\"?>\n";
            output << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
                      "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
            output << "  <UnstructuredGrid>\n";
            output << "    <Piece NumberOfPoints=\"" << s.cfg.npoin
                   << "\" NumberOfCells=\"" << total_cells << "\">\n";

            output << "      <Points>\n";
            output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                      "format=\"appended\" offset=\""
                   << block_offsets[points_block] << "\"/>\n";
            output << "      </Points>\n";

            output << "      <Cells>\n";
            output << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
                      "format=\"appended\" offset=\""
                   << block_offsets[connectivity_block] << "\"/>\n";
            output << "        <DataArray type=\"Int64\" Name=\"offsets\" "
                      "format=\"appended\" offset=\""
                   << block_offsets[offsets_block] << "\"/>\n";
            output << "        <DataArray type=\"UInt8\" Name=\"types\" "
                      "format=\"appended\" offset=\""
                   << block_offsets[types_block] << "\"/>\n";
            output << "      </Cells>\n";

            output << "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n";
            output << "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[pressure_block] << "\"/>\n";
            output << "        <DataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[temperature_block] << "\"/>\n";
            output << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
                   << block_offsets[velocity_block] << "\"/>\n";
            output << "        <DataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[velocity_magnitude_block] << "\"/>\n";
            output << "        <DataArray type=\"Int64\" Name=\"global_node_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[global_node_id_block] << "\"/>\n";
            output << "        <DataArray type=\"Int32\" Name=\"owner_rank\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[owner_rank_block] << "\"/>\n";
            output << "        <DataArray type=\"UInt8\" Name=\"is_owned\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[is_owned_block] << "\"/>\n";
            output << "        <DataArray type=\"UInt8\" Name=\"vtkGhostType\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[ghost_type_block] << "\"/>\n";
            output << "        <DataArray type=\"Int32\" Name=\"node_domain_kind\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[node_domain_kind_block] << "\"/>\n";
            output << "        <DataArray type=\"Int32\" Name=\"velocity_bc_type\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[velocity_bc_type_block] << "\"/>\n";
            output << "        <DataArray type=\"UInt8\" Name=\"pressure_fixed\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[pressure_fixed_block] << "\"/>\n";
            output << "      </PointData>\n";

            output << "      <CellData Scalars=\"material_id\">\n";
            output << "        <DataArray type=\"UInt8\" Name=\"cell_kind\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[cell_kind_block] << "\"/>\n";
            output << "        <DataArray type=\"Int64\" Name=\"global_element_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[global_element_id_block] << "\"/>\n";
            output << "        <DataArray type=\"Int32\" Name=\"material_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[material_id_block] << "\"/>\n";
            output << "        <DataArray type=\"Int32\" Name=\"bc_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[bc_id_block] << "\"/>\n";
            output << "        <DataArray type=\"Int64\" Name=\"parent_global_element\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                   << block_offsets[parent_global_element_block] << "\"/>\n";
            output << "      </CellData>\n";
            output << "    </Piece>\n";
            output << "  </UnstructuredGrid>\n";
            output << "  <AppendedData encoding=\"raw\">\n_";

            const auto write_block =
                [&](const VtuBlock block, const auto& fill_block)
                {
                    AppendedBlockWriter writer(
                        output,
                        block_sizes[static_cast<std::size_t>(block)]);
                    fill_block(writer);
                    writer.finish();
                };

            write_block(
                points_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        for (Int component = 1; component <= 3; ++component)
                        {
                            const double value = static_cast<double>(
                                safe_value(s.coord(component, ip)));
                            writer.append(value);
                        }
                    }
                });

            write_block(
                connectivity_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        for (Int node = 1; node <= 4; ++node)
                        {
                            writer.append(vtk_local_index(s.intma(node, ie)));
                        }
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        for (Int node = 1; node <= 3; ++node)
                        {
                            writer.append(vtk_local_index(s.iside(node, ib)));
                        }
                    }
                });

            write_block(
                offsets_block,
                [&](AppendedBlockWriter& writer)
                {
                    std::int64_t offset = 0;

                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        offset += 4;
                        writer.append(offset);
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        offset += 3;
                        writer.append(offset);
                    }
                });

            write_block(
                types_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(static_cast<std::uint8_t>(10));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        writer.append(static_cast<std::uint8_t>(5));
                    }
                });

            write_block(
                pressure_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        const double value = touches_fluid(s, ip)
                            ? static_cast<double>(safe_value(s.pres(ip)))
                            : 0.0;
                        writer.append(value);
                    }
                });

            write_block(
                temperature_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(static_cast<double>(
                            safe_value(s.temperature(ip))));
                    }
                });

            write_block(
                velocity_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        for (Int component = 1; component <= 3; ++component)
                        {
                            const double value = velocity_active(s, ip)
                                ? static_cast<double>(
                                    safe_value(s.unkno(component, ip)))
                                : 0.0;
                            writer.append(value);
                        }
                    }
                });

            write_block(
                velocity_magnitude_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        const double value = velocity_active(s, ip)
                            ? static_cast<double>(safe_value(s.velocity(ip)))
                            : 0.0;
                        writer.append(value);
                    }
                });

            write_block(
                global_node_id_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(static_cast<std::int64_t>(
                            global_node_id(s, ip)));
                    }
                });

            write_block(
                owner_rank_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(vtk_int32(
                            s.node_owner_rank[static_cast<std::size_t>(ip)],
                            "owner rank"));
                    }
                });

            write_block(
                is_owned_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(static_cast<std::uint8_t>(
                            s.node_owner_rank[static_cast<std::size_t>(ip)] ==
                                    s.mpi_rank
                                ? 1
                                : 0));
                    }
                });

            write_block(
                ghost_type_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(static_cast<std::uint8_t>(
                            s.node_owner_rank[static_cast<std::size_t>(ip)] ==
                                    s.mpi_rank
                                ? 0
                                : 1));
                    }
                });

            write_block(
                node_domain_kind_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(vtk_int32(
                            s.node_material_mask(ip),
                            "node domain kind"));
                    }
                });

            write_block(
                velocity_bc_type_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(vtk_int32(
                            s.node_velocity_bc_type(ip),
                            "velocity boundary type"));
                    }
                });

            write_block(
                pressure_fixed_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        writer.append(static_cast<std::uint8_t>(
                            s.node_pressure_fixed(ip) != 0 ? 1 : 0));
                    }
                });

            write_block(
                cell_kind_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(static_cast<std::uint8_t>(1));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        writer.append(static_cast<std::uint8_t>(2));
                    }
                });

            write_block(
                global_element_id_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(static_cast<std::int64_t>(
                            global_element_id(s, ie)));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        writer.append(static_cast<std::int64_t>(-1));
                    }
                });

            write_block(
                material_id_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(vtk_int32(
                            s.mat_elem(ie),
                            "material ID"));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        writer.append(static_cast<std::int32_t>(-1));
                    }
                });

            write_block(
                bc_id_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(static_cast<std::int32_t>(0));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        writer.append(vtk_int32(
                            s.iside(s.cfg.bsid, ib),
                            "boundary ID"));
                    }
                });

            write_block(
                parent_global_element_block,
                [&](AppendedBlockWriter& writer)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        writer.append(static_cast<std::int64_t>(
                            global_element_id(s, ie)));
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        const Int parent = s.iside(s.cfg.nsidpe, ib);
                        writer.append(static_cast<std::int64_t>(
                            global_element_id(s, parent)));
                    }
                });

            output << "\n  </AppendedData>\n";
            output << "</VTKFile>\n";
            output.flush();

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost failed while finalising a binary VTU piece");
            }
        }

'''

TIME_OLD = """        s.output_iterations.push_back(iteration);\n        s.output_times.push_back(s.cfg.rtime);\n"""

TIME_NEW = """        s.output_iterations.push_back(iteration);\n\n        const Real output_time = s.cfg.transient_on > 0\n            ? s.cfg.rtime\n            : static_cast<Real>(iteration);\n\n        s.output_times.push_back(output_time);\n"""


def replace_once(text, old, new, description):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(
            "{}: expected one match, found {}".format(description, count)
        )
    return text.replace(old, new, 1)


def main():
    if not SOURCE.is_file():
        raise FileNotFoundError("Missing source file: {}".format(SOURCE))

    text = SOURCE.read_text(encoding="utf-8")

    if "AppendedBlockWriter" in text and "header_type" in text:
        print("DD-5B1 binary VTU patch is already applied")
        return 0

    text = replace_once(text, INCLUDE_OLD, INCLUDE_NEW, "include block")

    start = text.find(START_MARKER)
    end = text.find(END_MARKER, start + 1)

    if start < 0 or end < 0 or end <= start:
        raise RuntimeError("Cannot locate the existing write_piece implementation")

    text = text[:start] + REPLACEMENT + text[end:]
    text = replace_once(text, TIME_OLD, TIME_NEW, "steady PVD time update")

    temporary = SOURCE.with_name(SOURCE.name + ".dd5b.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(SOURCE)

    print("DD-5B1 appended raw binary VTU patch: APPLIED")
    print("  modified: {}".format(SOURCE))
    print("  numerical solver loop: unchanged")
    print("  VTU encoding: appended raw binary, UInt64 headers")
    print("  steady PVD time: iteration number")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
