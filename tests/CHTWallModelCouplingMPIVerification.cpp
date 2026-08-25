#include "cbs/turbulence/CHTWallModelCoupling.hpp"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <exception>

#include <mpi.h>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;
    using cbs::turbulence::CHTWallModelCoupling;

    void set_face_map(CBSStateSI& s)
    {
        s.ippn1.resize(4, 3);
        s.ippn1(1, 1) = 2; s.ippn1(1, 2) = 3; s.ippn1(1, 3) = 4;
        s.ippn1(2, 1) = 1; s.ippn1(2, 2) = 4; s.ippn1(2, 3) = 3;
        s.ippn1(3, 1) = 1; s.ippn1(3, 2) = 2; s.ippn1(3, 3) = 4;
        s.ippn1(4, 1) = 1; s.ippn1(4, 2) = 3; s.ippn1(4, 3) = 2;
    }

    CBSStateSI make_rank_state(const Int rank)
    {
        CBSStateSI s;
        s.cfg.npoin = 4;
        s.cfg.nelem = 1;
        s.cfg.nboun = 0;
        s.cfg.temp_calc = 1;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.turbulent_thermal_diffusivity_on = 1;
        s.cfg.dimensional_mode = 1;
        s.cfg.material_properties_enabled = 1;
        s.cfg.sa_prandtl_t = 0.90;

        s.mpi_enabled = true;
        s.mpi_rank = rank;
        s.mpi_size = 2;

        set_face_map(s);

        s.intma.resize(4, 1);
        s.intma(1, 1) = 1;
        s.intma(2, 1) = 2;
        s.intma(3, 1) = 3;
        s.intma(4, 1) = 4;

        s.mat_elem.resize(1);
        s.mat_elem(1) = (rank == 0) ? 0 : 1;

        s.node_material_mask.resize(4);
        s.node_material_mask(1) = 3;
        s.node_material_mask(2) = 3;
        s.node_material_mask(3) = 3;
        s.node_material_mask(4) =
            (rank == 0)
                ? CBSStateSI::node_touches_fluid
                : CBSStateSI::node_touches_solid;

        // The face nodes are the same global nodes on both ranks.  The opposite
        // vertex differs, exactly as in a rank cut through a conformal interface.
        s.local_to_global_node.assign(5, 0);
        s.local_to_global_node[1] = 101;
        s.local_to_global_node[2] = 102;
        s.local_to_global_node[3] = 103;
        s.local_to_global_node[4] = (rank == 0) ? 104 : 105;

        s.detJ.resize(1);
        s.detJ(1) = 1.0e-3;

        s.annxf.resize(4, 4, 1);
        s.annxf.fill(0.0);
        s.annxf(3, 4, 1) = (rank == 0) ? -0.5 : 0.5;
        s.annxf(4, 4, 1) = 0.5;

        s.fedge.resize(4, 1);
        s.fedge.fill(0);

        s.dNkdx.resize(12);
        s.dNkdx.fill(0.0);
        const Real gradients[3][4] =
        {
            {-1.0, 1.0, 0.0, 0.0},
            {-1.0, 0.0, 1.0, 0.0},
            {-1000.0, 0.0, 0.0, 1000.0}
        };
        for (Int dim = 1; dim <= 3; ++dim)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                s.dNkdx((dim - 1) * 4 + a) = gradients[dim - 1][a - 1];
            }
        }

        s.delte.resize(1);
        s.delte(1) = 1.0e-7;

        s.rho_e.resize(1);
        s.cp_e.resize(1);
        s.rho_cp_e.resize(1);
        s.mu_e.resize(1);
        s.k_e.resize(1);
        s.nu_t_e.resize(1);
        s.k_eff_e.resize(1);

        if (rank == 0)
        {
            s.rho_e(1) = 6.7;
            s.cp_e(1) = 5200.0;
            s.rho_cp_e(1) = s.rho_e(1) * s.cp_e(1);
            s.mu_e(1) = 3.1e-5;
            s.k_e(1) = 0.24;
            s.nu_t_e(1) = 1.0e-5;
            s.k_eff_e(1) = s.k_e(1)
                + s.rho_cp_e(1) * s.nu_t_e(1) / s.cfg.sa_prandtl_t;
        }
        else
        {
            s.rho_e(1) = 7650.0;
            s.cp_e(1) = 560.0;
            s.rho_cp_e(1) = s.rho_e(1) * s.cp_e(1);
            s.mu_e(1) = 0.0;
            s.k_e(1) = 28.3;
            s.nu_t_e(1) = 0.0;
            s.k_eff_e(1) = 28.3;
        }

        s.unkn1.resize(3, 4);
        s.unkno.resize(3, 4);
        s.unkn1.fill(0.0);
        s.unkno.fill(0.0);
        if (rank == 0)
        {
            s.unkn1(1, 4) = 10.0;
            s.unkno(1, 4) = 10.0;
        }

        s.rhs.resize(3, 4);
        s.rhs.fill(0.0);

        s.temperature1.resize(4);
        s.temperature.resize(4);
        s.rhs1.resize(4);
        s.temperature1.fill(600.0);
        s.temperature.fill(600.0);
        s.rhs1.fill(0.0);

        s.node_velocity_bc_type.resize(4);
        s.node_velocity_bc_type.fill(CBSStateSI::velocity_bc_free);
        s.node_velocity_bc_type(1) = CBSStateSI::velocity_bc_noslip;
        s.node_velocity_bc_type(2) = CBSStateSI::velocity_bc_noslip;
        s.node_velocity_bc_type(3) = CBSStateSI::velocity_bc_noslip;

        return s;
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2)
    {
        if (rank == 0)
        {
            std::printf("FAIL: test requires exactly two MPI ranks\n");
        }
        MPI_Finalize();
        return 1;
    }

#if defined(_WIN32)
    _putenv_s("CBS3D_SA_WALL_TREATMENT", "");
    _putenv_s("CBS3D_CHT_WALL_TREATMENT", "1");
#else
    unsetenv("CBS3D_SA_WALL_TREATMENT");
    setenv("CBS3D_CHT_WALL_TREATMENT", "1", 1);
#endif

    int local_failure = 0;

    try
    {
        CBSStateSI s = make_rank_state(rank);

        const long long global_faces =
            CHTWallModelCoupling::globalWallFaceCount(s);

        if (global_faces != 1)
        {
            local_failure = 1;
        }

        const auto momentum =
            CHTWallModelCoupling::addMomentumWallFlux(s);

        if (rank == 0)
        {
            if (momentum.local_faces != 1 ||
                !(momentum.modeled_surface_load[0] < 0.0))
            {
                local_failure = 1;
            }
        }
        else if (momentum.local_faces != 0)
        {
            local_failure = 1;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("rank %d FAIL: %s\n", rank, error.what());
        local_failure = 1;
    }

    int global_failure = 0;
    MPI_Allreduce(
        &local_failure,
        &global_failure,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD);

    if (rank == 0 && global_failure == 0)
    {
        std::printf("PASS: two-rank conformal CHT wall reconstruction\n");
    }

    MPI_Finalize();
    return global_failure;
}
