//=============================================================================
// CBS3D++_SI
// Depth-dependent volumetric heat source (nuclear heating)
//=============================================================================

#include "cbs/preprocess/NuclearHeating.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    namespace
    {
        //=====================================================================
        // Length of the configured first-wall normal.
        //=====================================================================
        Real normal_length(const RunConfig& cfg)
        {
            return std::sqrt(
                cfg.nuclear_normal[0] * cfg.nuclear_normal[0] +
                cfg.nuclear_normal[1] * cfg.nuclear_normal[1] +
                cfg.nuclear_normal[2] * cfg.nuclear_normal[2]);
        }

        //=====================================================================
        // Index of the profile for one material, or -1 when none is given.
        //=====================================================================
        int profile_index(const RunConfig& cfg, const Int material_id)
        {
            for (std::size_t i = 0; i < cfg.nuclear_profiles.size(); ++i)
            {
                if (cfg.nuclear_profiles[i].material_id == material_id)
                {
                    return static_cast<int>(i);
                }
            }

            return -1;
        }

        //=====================================================================
        // Checks the profile list once, before any element is touched.
        //=====================================================================
        void validate_profiles(const RunConfig& cfg)
        {
            if (cfg.nuclear_profiles.empty())
            {
                throw std::runtime_error(
                    "NuclearHeating - nuclear_heating_on = 1 but no profile is listed");
            }

            const Real length = normal_length(cfg);

            if (!(length > 0.0) || !std::isfinite(length))
            {
                throw std::runtime_error(
                    "NuclearHeating - the first-wall normal (nx ny nz) must be non-zero");
            }

            for (std::size_t i = 0; i < cfg.nuclear_profiles.size(); ++i)
            {
                const NuclearHeatingProfile& p = cfg.nuclear_profiles[i];

                if (p.material_id < 0)
                {
                    throw std::runtime_error(
                        "NuclearHeating - profile material ID must be >= 0");
                }

                if (!std::isfinite(p.amplitude) || p.amplitude < 0.0)
                {
                    throw std::runtime_error(
                        "NuclearHeating - profile amplitude A must be finite and >= 0 "
                        "(material " + std::to_string(p.material_id) + ")");
                }

                if (!std::isfinite(p.decay_length) || !(p.decay_length > 0.0))
                {
                    throw std::runtime_error(
                        "NuclearHeating - profile decay length lambda must be > 0 "
                        "(material " + std::to_string(p.material_id) + ")");
                }

                for (std::size_t j = 0; j < i; ++j)
                {
                    if (cfg.nuclear_profiles[j].material_id == p.material_id)
                    {
                        throw std::runtime_error(
                            "NuclearHeating - material " + std::to_string(p.material_id)
                            + " has two profiles");
                    }
                }
            }
        }

        //=====================================================================
        // Effective element source, identical to the Step-4 assembly rule:
        //   fluid: Qvol_e;  solid: source_solid if non-zero, else Qvol_e.
        //=====================================================================
        Real effective_source(const CBSStateSI& s, const Int ie)
        {
            if (s.mat_elem(ie) != 0 && s.cfg.source_solid != 0.0)
            {
                return s.cfg.source_solid;
            }

            return s.Qvol_e(ie);
        }

        //=====================================================================
        // Sums an array over all ranks (no-op in a serial run).
        //=====================================================================
        void global_sum(const CBSStateSI& s, std::vector<double>& values)
        {
#ifdef CBS3D_USE_MPI
            if (s.mpi_enabled && !values.empty())
            {
                std::vector<double> total(values.size(), 0.0);

                if (MPI_Allreduce(
                        values.data(),
                        total.data(),
                        static_cast<int>(values.size()),
                        MPI_DOUBLE,
                        MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "NuclearHeating - MPI_Allreduce of the heat input failed");
                }

                values.swap(total);
            }
#else
            (void)s;
            (void)values;
#endif
        }
    }


    //=========================================================================
    // Q = A exp(-d / lambda), d = max(0, (x - x0) . n / |n|)
    //=========================================================================
    Real NuclearHeating::profileValue(
        const RunConfig& cfg,
        const NuclearHeatingProfile& profile,
        const Real x,
        const Real y,
        const Real z)
    {
        const Real length = normal_length(cfg);

        Real depth =
            ((x - cfg.nuclear_origin[0]) * cfg.nuclear_normal[0] +
             (y - cfg.nuclear_origin[1]) * cfg.nuclear_normal[1] +
             (z - cfg.nuclear_origin[2]) * cfg.nuclear_normal[2]) / length;

        if (depth < 0.0)
        {
            depth = 0.0;
        }

        return profile.amplitude * std::exp(-depth / profile.decay_length);
    }


    //=========================================================================
    // Writes the profile value at every element centroid into Qvol_e.
    //=========================================================================
    void NuclearHeating::applyProfiles(CBSStateSI& s)
    {
        if (s.cfg.nuclear_heating_on < 1)
        {
            return;
        }

        if (s.cfg.material_properties_enabled < 1)
        {
            throw std::runtime_error(
                "NuclearHeating - the nuclear-heating block needs "
                "material_properties_enabled = 1 (.material and .matprop files)");
        }

        validate_profiles(s.cfg);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Int material = s.mat_elem(ie);
            const int ip = profile_index(s.cfg, material);

            if (ip < 0)
            {
                continue;
            }

            if (s.Qvol_e(ie) != 0.0)
            {
                throw std::runtime_error(
                    "NuclearHeating - material " + std::to_string(material)
                    + " has both a constant .matprop Qvol and a nuclear-heating "
                    "profile; set its .matprop Qvol to 0");
            }

            if (material != 0 && s.cfg.source_solid != 0.0)
            {
                throw std::runtime_error(
                    "NuclearHeating - .par source_solid is non-zero and solid material "
                    + std::to_string(material)
                    + " also has a nuclear-heating profile; set source_solid to 0");
            }

            Real xc = 0.0;
            Real yc = 0.0;
            Real zc = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int node = s.intma(a, ie);
                xc += s.coord(1, node);
                yc += s.coord(2, node);
                zc += s.coord(3, node);
            }

            const Real inv = 1.0 / static_cast<Real>(s.cfg.nep);

            s.Qvol_e(ie) = profileValue(
                s.cfg,
                s.cfg.nuclear_profiles[static_cast<std::size_t>(ip)],
                xc * inv,
                yc * inv,
                zc * inv);
        }
    }


    //=========================================================================
    // Global heat input: volumetric (per profile material and the rest) and
    // the BC 532 surface flux. Printed by rank 0.
    //
    // Layout of the reduced array:
    //     [3 i + 0] Q   of profile i   [W]
    //     [3 i + 1] V   of profile i   [m^3]
    //     [3 i + 2] elements of profile i
    //     [3 n + 0] Q   of all other elements (constant sources)
    //     [3 n + 1] area of BC 532 faces  [m^2]
    //=========================================================================
    void NuclearHeating::reportHeatInput(const CBSStateSI& s)
    {
        if (s.cfg.temp_calc < 1)
        {
            return;
        }

        const std::size_t n =
            s.cfg.nuclear_heating_on > 0 ? s.cfg.nuclear_profiles.size() : 0;

        std::vector<double> sums(3 * n + 2, 0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            // detJ = 6 V for a P1 tetrahedron.
            const Real volume = s.detJ(ie) / 6.0;
            const Real q = effective_source(s, ie);
            const int ip = n > 0 ? profile_index(s.cfg, s.mat_elem(ie)) : -1;

            if (ip >= 0)
            {
                sums[3 * static_cast<std::size_t>(ip) + 0] += q * volume;
                sums[3 * static_cast<std::size_t>(ip) + 1] += volume;
                sums[3 * static_cast<std::size_t>(ip) + 2] += 1.0;
            }
            else
            {
                sums[3 * n + 0] += q * volume;
            }
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            if (s.iside(s.cfg.bsid, ib) == s.cfg.bc_noslip_heatflux_wall)
            {
                sums[3 * n + 1] += s.face_norm(4, ib);
            }
        }

        global_sum(s, sums);

        for (std::size_t i = 0; i < n; ++i)
        {
            if (sums[3 * i + 2] < 0.5)
            {
                throw std::runtime_error(
                    "NuclearHeating - profile material "
                    + std::to_string(s.cfg.nuclear_profiles[i].material_id)
                    + " has no elements in the mesh (check the .matprop IDs)");
            }
        }

        if (s.mpi_rank != 0)
        {
            return;
        }

        const double q_surface = s.cfg.heat_flux_bc * sums[3 * n + 1];
        double q_volume = sums[3 * n + 0];

        // Written with C stdio, not std::cout: the MPI production loop
        // silences std::cout during set-up, and the energy-balance check
        // (tools/energy_balance.py) reads HEAT_INPUT_TOTAL_W from the log.
        std::ostringstream out;
        out << std::scientific << std::setprecision(6);
        out << "  HEAT INPUT (whole domain)\n";

        for (std::size_t i = 0; i < n; ++i)
        {
            const NuclearHeatingProfile& p = s.cfg.nuclear_profiles[i];
            q_volume += sums[3 * i + 0];

            out
                << "    nuclear heating, material " << p.material_id
                << " : A = " << p.amplitude << " W/m3, lambda = "
                << p.decay_length << " m, volume = " << sums[3 * i + 1]
                << " m3, Q = " << sums[3 * i + 0] << " W\n";
        }

        if (sums[3 * n + 0] != 0.0)
        {
            out
                << "    constant volumetric sources (.matprop Qvol / source_solid) : Q = "
                << sums[3 * n + 0] << " W\n";
        }

        out
            << "    BC 532 surface flux : q = " << s.cfg.heat_flux_bc
            << " W/m2, area = " << sums[3 * n + 1]
            << " m2, Q = " << q_surface << " W\n"
            << "    HEAT_INPUT_TOTAL_W = " << (q_volume + q_surface) << "\n";

        std::fputs(out.str().c_str(), stdout);
        std::fflush(stdout);
    }
}
