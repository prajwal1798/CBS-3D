//=============================================================================
// CBS3D++_SI
//
// Standalone boundedness test for the production Step-4 thermal transport
// operator on a small structured tetrahedral box.
//
// The test advances a bounded 300--400 K top-hat temperature profile with a
// prescribed constant velocity.  No cooling source is present.  Therefore the
// continuous advection-diffusion problem must satisfy:
//
//     300 K <= T(x,t) <= 400 K.
//
// EnergyAssembly::assembleStep4Rhs() is called directly at every step.  The
// program reports bound violations and excess-energy conservation for two time
// steps at the same final time.  It is a diagnostic characterisation executable
// and deliberately does not alter the production formulation.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    constexpr Real reference_temperature = 300.0;
    constexpr Real maximum_initial_temperature = 400.0;
    constexpr Real bound_tolerance = 1.0e-10;

    constexpr Real domain_length_x = 2.0e-2;
    constexpr Real spacing = 4.0e-4;
    constexpr Int nodes_x = 51;
    constexpr Int nodes_y = 3;
    constexpr Int nodes_z = 3;

    constexpr Real velocity_x = 5.0e-2;
    constexpr Real density = 997.0;
    constexpr Real heat_capacity = 4182.0;
    constexpr Real conductivity = 0.606;
    constexpr Real rho_cp = density * heat_capacity;
    constexpr Real thermal_diffusivity = conductivity / rho_cp;

    constexpr Real pulse_start = 6.0e-3;
    constexpr Real pulse_end = 8.0e-3;
    constexpr Real final_time = 1.0e-1;

    struct Mesh
    {
        std::vector<std::array<Real, 3>> coordinates;
        std::vector<std::array<Int, 4>> tetrahedra;
    };

    struct Metrics
    {
        Real minimum = std::numeric_limits<Real>::max();
        Real maximum = std::numeric_limits<Real>::lowest();
        Int below_count = 0;
        Int above_count = 0;
    };

    Int node_id(const Int i, const Int j, const Int k)
    {
        return 1 + i + nodes_x * (j + nodes_y * k);
    }

    Real determinant3(const std::array<std::array<Real, 3>, 3>& matrix)
    {
        return
            matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
          - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
          + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    }

    std::array<std::array<Real, 3>, 3> inverse3(
        const std::array<std::array<Real, 3>, 3>& matrix,
        const Real determinant)
    {
        if (std::abs(determinant) <= 1.0e-30)
        {
            throw std::runtime_error("Step4BoundedTransport - singular tetrahedron");
        }

        std::array<std::array<Real, 3>, 3> inverse{};

        inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant;
        inverse[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant;
        inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant;

        inverse[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant;
        inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant;
        inverse[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant;

        inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant;
        inverse[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant;
        inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant;

        return inverse;
    }

    std::array<std::array<Real, 3>, 3> jacobian(
        const Mesh& mesh,
        const std::array<Int, 4>& tetrahedron)
    {
        const auto& x0 = mesh.coordinates[static_cast<std::size_t>(tetrahedron[0] - 1)];
        const auto& x1 = mesh.coordinates[static_cast<std::size_t>(tetrahedron[1] - 1)];
        const auto& x2 = mesh.coordinates[static_cast<std::size_t>(tetrahedron[2] - 1)];
        const auto& x3 = mesh.coordinates[static_cast<std::size_t>(tetrahedron[3] - 1)];

        std::array<std::array<Real, 3>, 3> matrix{};
        for (Int dim = 0; dim < 3; ++dim)
        {
            matrix[static_cast<std::size_t>(dim)][0] = x1[static_cast<std::size_t>(dim)] - x0[static_cast<std::size_t>(dim)];
            matrix[static_cast<std::size_t>(dim)][1] = x2[static_cast<std::size_t>(dim)] - x0[static_cast<std::size_t>(dim)];
            matrix[static_cast<std::size_t>(dim)][2] = x3[static_cast<std::size_t>(dim)] - x0[static_cast<std::size_t>(dim)];
        }
        return matrix;
    }

    void append_oriented_tetrahedron(
        Mesh& mesh,
        std::array<Int, 4> tetrahedron)
    {
        Real determinant = determinant3(jacobian(mesh, tetrahedron));
        if (determinant < 0.0)
        {
            std::swap(tetrahedron[1], tetrahedron[2]);
            determinant = determinant3(jacobian(mesh, tetrahedron));
        }

        if (determinant <= 0.0)
        {
            throw std::runtime_error("Step4BoundedTransport - invalid tetrahedral orientation");
        }

        mesh.tetrahedra.push_back(tetrahedron);
    }

    Mesh make_mesh()
    {
        Mesh mesh;
        mesh.coordinates.resize(static_cast<std::size_t>(nodes_x * nodes_y * nodes_z));

        for (Int k = 0; k < nodes_z; ++k)
        {
            for (Int j = 0; j < nodes_y; ++j)
            {
                for (Int i = 0; i < nodes_x; ++i)
                {
                    const Int id = node_id(i, j, k);
                    mesh.coordinates[static_cast<std::size_t>(id - 1)] =
                    {
                        spacing * static_cast<Real>(i),
                        spacing * static_cast<Real>(j),
                        spacing * static_cast<Real>(k)
                    };
                }
            }
        }

        for (Int k = 0; k < nodes_z - 1; ++k)
        {
            for (Int j = 0; j < nodes_y - 1; ++j)
            {
                for (Int i = 0; i < nodes_x - 1; ++i)
                {
                    const Int v000 = node_id(i,     j,     k);
                    const Int v100 = node_id(i + 1, j,     k);
                    const Int v010 = node_id(i,     j + 1, k);
                    const Int v110 = node_id(i + 1, j + 1, k);
                    const Int v001 = node_id(i,     j,     k + 1);
                    const Int v101 = node_id(i + 1, j,     k + 1);
                    const Int v011 = node_id(i,     j + 1, k + 1);
                    const Int v111 = node_id(i + 1, j + 1, k + 1);

                    append_oriented_tetrahedron(mesh, {v000, v100, v110, v111});
                    append_oriented_tetrahedron(mesh, {v000, v110, v010, v111});
                    append_oriented_tetrahedron(mesh, {v000, v010, v011, v111});
                    append_oriented_tetrahedron(mesh, {v000, v011, v001, v111});
                    append_oriented_tetrahedron(mesh, {v000, v001, v101, v111});
                    append_oriented_tetrahedron(mesh, {v000, v101, v100, v111});
                }
            }
        }

        return mesh;
    }

    Int gradient_index(
        const CBSStateSI& state,
        const Int element,
        const Int dimension,
        const Int local_node)
    {
        return (element - 1) * state.cfg.ndim * state.cfg.nep
             + (dimension - 1) * state.cfg.nep
             + local_node;
    }

    void populate_geometry(CBSStateSI& state, const Mesh& mesh)
    {
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            const auto& coordinate = mesh.coordinates[static_cast<std::size_t>(node - 1)];
            for (Int dim = 1; dim <= 3; ++dim)
            {
                state.coord(dim, node) = coordinate[static_cast<std::size_t>(dim - 1)];
            }
        }

        for (Int element = 1; element <= state.cfg.nelem; ++element)
        {
            const auto& tetrahedron = mesh.tetrahedra[static_cast<std::size_t>(element - 1)];
            for (Int local_node = 1; local_node <= 4; ++local_node)
            {
                state.intma(local_node, element) = tetrahedron[static_cast<std::size_t>(local_node - 1)];
            }

            const auto matrix = jacobian(mesh, tetrahedron);
            const Real determinant = determinant3(matrix);
            const auto inverse = inverse3(matrix, determinant);

            state.detJ(element) = determinant;

            std::array<std::array<Real, 3>, 4> gradients{};
            gradients[1] = {inverse[0][0], inverse[0][1], inverse[0][2]};
            gradients[2] = {inverse[1][0], inverse[1][1], inverse[1][2]};
            gradients[3] = {inverse[2][0], inverse[2][1], inverse[2][2]};

            for (Int dim = 0; dim < 3; ++dim)
            {
                gradients[0][static_cast<std::size_t>(dim)] =
                    -gradients[1][static_cast<std::size_t>(dim)]
                    -gradients[2][static_cast<std::size_t>(dim)]
                    -gradients[3][static_cast<std::size_t>(dim)];
            }

            for (Int local_node = 1; local_node <= 4; ++local_node)
            {
                for (Int dim = 1; dim <= 3; ++dim)
                {
                    state.dNkdx(gradient_index(state, element, dim, local_node)) =
                        gradients[static_cast<std::size_t>(local_node - 1)]
                                 [static_cast<std::size_t>(dim - 1)];
                }
            }
        }
    }

    CBSStateSI make_state(const Mesh& mesh, const Real dt)
    {
        CBSStateSI state;
        state.set_problem_sizes(
            static_cast<Int>(mesh.tetrahedra.size()),
            static_cast<Int>(mesh.coordinates.size()),
            0,
            0);

        state.cfg.temp_calc = 1;
        state.cfg.heat_flux_bc = 0.0;
        state.cfg.turbulence_on = 0;
        state.cfg.turbulent_thermal_diffusivity_on = 0;

        populate_geometry(state, mesh);

        for (Int element = 1; element <= state.cfg.nelem; ++element)
        {
            state.mat_elem(element) = 0;
            state.delte(element) = dt;
            state.rho_e(element) = density;
            state.cp_e(element) = heat_capacity;
            state.rho_cp_e(element) = rho_cp;
            state.k_e(element) = conductivity;
            state.Qvol_e(element) = 0.0;
        }

        state.unkno.fill(0.0);
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            state.unkno(1, node) = velocity_x;

            const Real x = state.coord(1, node);
            const Real temperature =
                (x >= pulse_start - 1.0e-14 && x <= pulse_end + 1.0e-14)
                ? maximum_initial_temperature
                : reference_temperature;

            state.temperature(node) = temperature;
            state.temperature1(node) = temperature;
        }

        return state;
    }

    std::vector<Real> assemble_lumped_capacity(
        CBSStateSI& state,
        const Real dt)
    {
        std::vector<Real> capacity(static_cast<std::size_t>(state.cfg.npoin + 1), 0.0);

        for (Int element = 1; element <= state.cfg.nelem; ++element)
        {
            const Real element_nodal_capacity =
                state.rho_cp_e(element) * state.detJ(element) / 24.0;

            for (Int local_node = 1; local_node <= state.cfg.nep; ++local_node)
            {
                const Int node = state.intma(local_node, element);
                capacity[static_cast<std::size_t>(node)] += element_nodal_capacity;
            }
        }

        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            const Real nodal_capacity = capacity[static_cast<std::size_t>(node)];
            if (nodal_capacity <= 0.0)
            {
                throw std::runtime_error("Step4BoundedTransport - zero nodal capacity");
            }
            state.elcoe2p(node) = dt / nodal_capacity;
        }

        return capacity;
    }

    void enforce_inlet_temperature(CBSStateSI& state)
    {
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            if (std::abs(state.coord(1, node)) <= 1.0e-14)
            {
                state.temperature(node) = reference_temperature;
            }
        }
    }

    Metrics measure(const CBSStateSI& state)
    {
        Metrics metrics;
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            const Real temperature = state.temperature(node);
            metrics.minimum = std::min(metrics.minimum, temperature);
            metrics.maximum = std::max(metrics.maximum, temperature);

            if (temperature < reference_temperature - bound_tolerance)
            {
                ++metrics.below_count;
            }
            if (temperature > maximum_initial_temperature + bound_tolerance)
            {
                ++metrics.above_count;
            }
        }
        return metrics;
    }

    Real excess_energy(
        const CBSStateSI& state,
        const std::vector<Real>& capacity)
    {
        Real energy = 0.0;
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            energy += capacity[static_cast<std::size_t>(node)]
                    * (state.temperature(node) - reference_temperature);
        }
        return energy;
    }

    void run_case(const Mesh& mesh, const Real dt)
    {
        const Int steps = static_cast<Int>(std::llround(final_time / dt));
        CBSStateSI state = make_state(mesh, dt);
        const std::vector<Real> capacity = assemble_lumped_capacity(state, dt);

        const Real initial_energy = excess_energy(state, capacity);
        Int first_violation_step = -1;

        for (Int step = 1; step <= steps; ++step)
        {
            for (Int node = 1; node <= state.cfg.npoin; ++node)
            {
                state.temperature1(node) = state.temperature(node);
            }

            cbs::EnergyAssembly::assembleStep4Rhs(state);

            for (Int node = 1; node <= state.cfg.npoin; ++node)
            {
                state.temperature(node) =
                    state.temperature1(node)
                    + state.rhs1(node) * state.elcoe2p(node);
            }

            enforce_inlet_temperature(state);

            if (first_violation_step < 0)
            {
                const Metrics current = measure(state);
                if (current.below_count > 0 || current.above_count > 0)
                {
                    first_violation_step = step;
                }
            }
        }

        const Metrics final_metrics = measure(state);
        const Real final_energy = excess_energy(state, capacity);
        const Real relative_energy_drift =
            (initial_energy != 0.0)
            ? (final_energy - initial_energy) / initial_energy
            : 0.0;

        const Real peclet = velocity_x * spacing / (2.0 * thermal_diffusivity);
        const Real courant = velocity_x * dt / spacing;

        std::cout << "\nCASE dt = " << std::scientific << dt << " s\n";
        std::cout << "  steps                         : " << steps << "\n";
        std::cout << "  final time                    : " << final_time << " s\n";
        std::cout << "  element Peclet U*h/(2*alpha) : " << peclet << "\n";
        std::cout << "  Courant U*dt/h                : " << courant << "\n";
        std::cout << "  Tmin                          : " << std::setprecision(12) << final_metrics.minimum << " K\n";
        std::cout << "  Tmax                          : " << final_metrics.maximum << " K\n";
        std::cout << "  nodes below 300 K             : " << final_metrics.below_count << "\n";
        std::cout << "  nodes above 400 K             : " << final_metrics.above_count << "\n";
        std::cout << "  first violation step          : " << first_violation_step << "\n";
        if (first_violation_step >= 0)
        {
            std::cout << "  first violation time          : "
                      << static_cast<Real>(first_violation_step) * dt << " s\n";
        }
        std::cout << "  initial excess energy         : " << initial_energy << " J\n";
        std::cout << "  final excess energy           : " << final_energy << " J\n";
        std::cout << "  relative excess-energy drift  : " << relative_energy_drift << "\n";
        std::cout << "  BOUND CHECK                   : "
                  << ((final_metrics.below_count == 0 && final_metrics.above_count == 0)
                      ? "PASS" : "FAIL")
                  << "\n";
    }
}

int main()
{
    const Mesh mesh = make_mesh();

    std::cout << "CBS3D++_SI Step-4 bounded thermal-transport validation\n";
    std::cout << "Production routine: EnergyAssembly::assembleStep4Rhs\n";
    std::cout << "Mesh nodes      : " << mesh.coordinates.size() << "\n";
    std::cout << "Mesh tetrahedra : " << mesh.tetrahedra.size() << "\n";
    std::cout << "Domain          : [0," << domain_length_x << "] x [0,"
              << spacing * static_cast<Real>(nodes_y - 1) << "] x [0,"
              << spacing * static_cast<Real>(nodes_z - 1) << "] m\n";
    std::cout << "Initial bounds  : [300,400] K\n";
    std::cout << "Prescribed U    : (0.05,0,0) m/s\n";

    run_case(mesh, 1.0e-4);
    run_case(mesh, 1.0e-5);

    std::cout << "\n============================================================\n";
    std::cout << "This diagnostic reports whether the current production Step-4\n";
    std::cout << "operator preserves the physical 300--400 K bounds.\n";
    return 0;
}
