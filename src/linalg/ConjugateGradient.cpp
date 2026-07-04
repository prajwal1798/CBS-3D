//=============================================================================
// CBS3D++_SI
//
// Serial Conjugate Gradient solver for the CBS pressure equation.
//
// The constrained linear system is:
//
//     A p = b
//
// The pressure unknown is active only at nodes connected to fluid elements.
// Prescribed-pressure nodes are fixed, and solid-only nodes are removed from
// all residuals and search directions.
//
// For preconditioned Conjugate Gradient:
//
//     r_0 = b - A p_0
//
//     M z_k = r_k
//
//     alpha_k = (r_k, z_k) / (d_k, A d_k)
//
//     p_(k+1) = p_k + alpha_k d_k
//
//     r_(k+1) = r_k - alpha_k A d_k
//
//     beta_k = (r_(k+1), z_(k+1)) / (r_k, z_k)
//
//     d_(k+1) = z_(k+1) + beta_k d_k
//
// With Jacobi preconditioning:
//
//     z_i = r_i / A_ii
//
// The matrix-vector product is evaluated through MatrixVectorCalc using the
// compact tetrahedral pressure operator.
//=============================================================================

#include "cbs/linalg/ConjugateGradient.hpp"

#include "cbs/linalg/MatrixVectorCalc.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        // Returns elapsed wall-clock time in seconds.
        double elapsed_seconds(
            const std::chrono::steady_clock::time_point start,
            const std::chrono::steady_clock::time_point finish)
        {
            return std::chrono::duration<double>(finish - start).count();
        }

        // Returns true when a node belongs to the fluid pressure space.
        bool is_pressure_active(
            const std::vector<char>& active,
            const Int ip)
        {
            return active[static_cast<std::size_t>(ip)] != 0;
        }

        // Returns true when pressure is prescribed at the node.
        bool is_fixed_pressure_node(
            const std::vector<char>& fixed,
            const Int ip)
        {
            return fixed[static_cast<std::size_t>(ip)] != 0;
        }

        // A free pressure node is fluid-connected and not constrained:
        //
        //     free = active AND NOT fixed
        bool is_free_pressure_node(
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Int ip)
        {
            return is_pressure_active(active, ip) &&
                   !is_fixed_pressure_node(fixed, ip);
        }

        // Builds the pressure-space mask from fluid-element connectivity.
        //
        // A node is active when it belongs to at least one element satisfying:
        //
        //     mat_elem(e) = 0
        //
        // Solid-only nodes are excluded from the pressure solve.
        void build_pressure_active_mask(
            const CBSStateSI& s,
            std::vector<char>& active)
        {
            active.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "ConjugateGradient::solvePressure - element connectivity node out of range");
                    }

                    active[static_cast<std::size_t>(ip)] = 1;
                }
            }

            Int active_count = 0;
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (is_pressure_active(active, ip))
                {
                    ++active_count;
                }
            }

            if (active_count == 0)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - no fluid-connected pressure nodes found");
            }
        }

        // Marks all prescribed-pressure and reference-pressure nodes.
        //
        // For every entry in bc_list:
        //
        //     fixed(node)        = 1
        //     fixed_values(node) = prescribed pressure
        void build_fixed_pressure_mask(
            const CBSStateSI& s,
            const std::vector<char>& active,
            std::vector<char>& fixed,
            Array1D<Real>& fixed_values)
        {
            fixed.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            fixed_values.resize(s.cfg.npoin);
            fixed_values.fill(0.0);

            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                const Int ip = s.bc_list(i);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "ConjugateGradient::solvePressure - fixed pressure node out of range");
                }

                if (!is_pressure_active(active, ip))
                {
                    throw std::runtime_error(
                        "ConjugateGradient::solvePressure - fixed pressure node is not fluid-connected: node "
                        + std::to_string(ip));
                }

                fixed[static_cast<std::size_t>(ip)] = 1;
                fixed_values(ip) = s.bc_values(i);
            }
        }

        // Applies pressure-space constraints to the solver arrays.
        //
        // Inactive solid-only node:
        //
        //     p_i = 0,  b_i = 0,  A_ii = 0
        //
        // Fixed pressure node:
        //
        //     p_i = p_bc,  b_i = p_bc,  A_ii = 1
        void impose_pressure_constraints(
            CBSStateSI& s,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Array1D<Real>& fixed_values)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (!is_pressure_active(active, ip))
                {
                    // Solid-only nodes do not belong to the pressure space.
                    // They are not pressure Dirichlet nodes; they are inactive.
                    s.pres(ip) = 0.0;
                    s.rhs1(ip) = 0.0;
                    s.pdiag(ip) = 0.0;
                    continue;
                }

                if (is_fixed_pressure_node(fixed, ip))
                {
                    // Physical/reference pressure rows are identity rows.
                    s.pres(ip) = fixed_values(ip);
                    s.rhs1(ip) = fixed_values(ip);
                    s.pdiag(ip) = 1.0;
                }
            }
        }

        // Removes inactive and prescribed-pressure entries from a CG
        // residual, preconditioned residual or search direction.
        void zero_inactive_or_fixed_nodes(
            Array1D<Real>& x,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Int npoin)
        {
            for (Int ip = 1; ip <= npoin; ++ip)
            {
                if (!is_free_pressure_node(active, fixed, ip))
                {
                    x(ip) = 0.0;
                }
            }
        }

        // Calculates the Euclidean inner product over free pressure nodes:
        //
        //     (a,b) = sum_i a_i b_i
        Real dot_owned_serial(
            const Array1D<Real>& a,
            const Array1D<Real>& b,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Int npoin)
        {
            Real sum = 0.0;

#pragma omp parallel for reduction(+:sum) schedule(static)
            for (Int ip = 1; ip <= npoin; ++ip)
            {
                if (is_free_pressure_node(active, fixed, ip))
                {
                    sum += a(ip) * b(ip);
                }
            }

            return sum;
        }

        // Calculates the Euclidean residual norm:
        //
        //     ||r||_2 = sqrt((r,r))
        Real l2_norm_serial(
            const Array1D<Real>& r,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Int npoin)
        {
            return std::sqrt(dot_owned_serial(r, r, active, fixed, npoin));
        }

        // Calculates the maximum absolute residual over free nodes:
        //
        //     ||r||_infinity = max_i |r_i|
        Real max_abs_serial(
            const Array1D<Real>& r,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Int npoin)
        {
            Real maxv = 0.0;

            for (Int ip = 1; ip <= npoin; ++ip)
            {
                if (is_free_pressure_node(active, fixed, ip))
                {
                    maxv = std::max(maxv, std::abs(r(ip)));
                }
            }

            return maxv;
        }

        // Applies the selected pressure preconditioner.
        //
        // No preconditioner:
        //
        //     z = r
        //
        // Jacobi preconditioner:
        //
        //     z_i = r_i / A_ii
        //
        // where A_ii is stored in pdiag.
        void apply_preconditioner(
            const CBSStateSI& s,
            const Array1D<Real>& r,
            Array1D<Real>& z,
            const std::vector<char>& active,
            const std::vector<char>& fixed)
        {
            z.fill(0.0);

            if (s.cfg.cg_preconditioner == 0)
            {
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (is_free_pressure_node(active, fixed, ip))
                    {
                        if (!std::isfinite(r(ip)))
                        {
                            throw std::runtime_error(
                                "ConjugateGradient::solvePressure - invalid pressure residual at node "
                                + std::to_string(ip)
                                + ", r="
                                + std::to_string(r(ip)));
                        }

                        z(ip) = r(ip);
                    }
                }

                return;
            }

            if (s.cfg.cg_preconditioner != 1)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - cg_preconditioner must be 0 or 1");
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (!is_free_pressure_node(active, fixed, ip))
                {
                    z(ip) = 0.0;
                    continue;
                }

                const Real diag = s.pdiag(ip);

                if (!std::isfinite(diag) || diag <= 1.0e-30)
                {
                    throw std::runtime_error(
                        "ConjugateGradient::solvePressure - invalid Jacobi pressure diagonal at active node "
                        + std::to_string(ip)
                        + ", pdiag="
                        + std::to_string(diag));
                }

                if (!std::isfinite(r(ip)))
                {
                    throw std::runtime_error(
                        "ConjugateGradient::solvePressure - invalid pressure residual at node "
                        + std::to_string(ip)
                        + ", r="
                        + std::to_string(r(ip)));
                }

                z(ip) = r(ip) / diag;
            }
        }

        // Collects diagnostic information when the preconditioned
        // residual inner product is invalid.
        std::string diagnose_initial_cg_state(
            const CBSStateSI& s,
            const Array1D<Real>& r,
            const Array1D<Real>& z,
            const std::vector<char>& active,
            const std::vector<char>& fixed,
            const Real rz_old)
        {
            Int bad_pdiag_node = 0;
            Int negative_pdiag_node = 0;
            Int bad_rhs_node = 0;
            Int bad_residual_node = 0;
            Int bad_z_node = 0;
            Int active_free_count = 0;

            Real min_pdiag = std::numeric_limits<Real>::max();
            Real max_pdiag = 0.0;
            Real max_abs_rhs = 0.0;
            Real max_abs_r = 0.0;
            Real max_abs_z = 0.0;

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (!is_free_pressure_node(active, fixed, ip))
                {
                    continue;
                }

                ++active_free_count;

                const Real diag = s.pdiag(ip);

                if (bad_pdiag_node == 0 &&
                    (!std::isfinite(diag) || std::abs(diag) <= 1.0e-300))
                {
                    bad_pdiag_node = ip;
                }

                if (negative_pdiag_node == 0 &&
                    std::isfinite(diag) &&
                    diag < 0.0)
                {
                    negative_pdiag_node = ip;
                }

                if (bad_rhs_node == 0 && !std::isfinite(s.rhs1(ip)))
                {
                    bad_rhs_node = ip;
                }

                if (bad_residual_node == 0 && !std::isfinite(r(ip)))
                {
                    bad_residual_node = ip;
                }

                if (bad_z_node == 0 && !std::isfinite(z(ip)))
                {
                    bad_z_node = ip;
                }

                if (std::isfinite(diag))
                {
                    min_pdiag = std::min(min_pdiag, diag);
                    max_pdiag = std::max(max_pdiag, diag);
                }

                if (std::isfinite(s.rhs1(ip)))
                {
                    max_abs_rhs = std::max(max_abs_rhs, std::abs(s.rhs1(ip)));
                }

                if (std::isfinite(r(ip)))
                {
                    max_abs_r = std::max(max_abs_r, std::abs(r(ip)));
                }

                if (std::isfinite(z(ip)))
                {
                    max_abs_z = std::max(max_abs_z, std::abs(z(ip)));
                }
            }

            if (active_free_count == 0)
            {
                min_pdiag = 0.0;
            }

            return
                "rz_old=" + std::to_string(rz_old) +
                ", active_free_count=" + std::to_string(active_free_count) +
                ", bad_pdiag_node=" + std::to_string(bad_pdiag_node) +
                ", negative_pdiag_node=" + std::to_string(negative_pdiag_node) +
                ", bad_rhs_node=" + std::to_string(bad_rhs_node) +
                ", bad_residual_node=" + std::to_string(bad_residual_node) +
                ", bad_z_node=" + std::to_string(bad_z_node) +
                ", min_pdiag=" + std::to_string(min_pdiag) +
                ", max_pdiag=" + std::to_string(max_pdiag) +
                ", max_abs_rhs=" + std::to_string(max_abs_rhs) +
                ", max_abs_r=" + std::to_string(max_abs_r) +
                ", max_abs_z=" + std::to_string(max_abs_z);
        }

        // Applies the selected pressure convergence test.
        //
        //     cg_conv_test = 1 : max_i |r_i| <= absToler
        //
        //     cg_conv_test = 2 : ||r_k||_2 / ||r_0||_2 <= relToler
        //
        //     cg_conv_test = 3 : either criterion is sufficient
        bool convergence_reached(
            const CBSStateSI& s,
            const Real rel_l2,
            const Real max_abs)
        {
            if (s.cfg.cg_conv_test == 1)
            {
                return max_abs <= s.cfg.absToler;
            }

            if (s.cfg.cg_conv_test == 2)
            {
                return rel_l2 <= s.cfg.relToler;
            }

            if (s.cfg.cg_conv_test == 3)
            {
                return max_abs <= s.cfg.absToler || rel_l2 <= s.cfg.relToler;
            }

            throw std::runtime_error(
                "ConjugateGradient::solvePressure - cg_conv_test must be 1, 2, or 3");
        }

        // Checks the Pressure CG iteration limit, tolerances,
        // preconditioner option and convergence-test option.
        void validate_cg_config(const CBSStateSI& s)
        {
            if (s.cfg.cg_max_iter < 1)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - cg_max_iter must be positive");
            }

            if (s.cfg.relToler <= 0.0 || s.cfg.absToler <= 0.0)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - relToler and absToler must be positive");
            }

            if (s.cfg.cg_preconditioner < 0 || s.cfg.cg_preconditioner > 1)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - cg_preconditioner must be 0 or 1");
            }

            if (s.cfg.cg_conv_test < 1 || s.cfg.cg_conv_test > 3)
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - cg_conv_test must be 1, 2, or 3");
            }
        }
    }

    //=========================================================================
    // Solves the constrained fluid-pressure system using preconditioned
    // Conjugate Gradient iteration.
    //
    // Linear system:
    //
    //     A p = b
    //
    // Initial residual:
    //
    //     r_0 = b - A p_0
    //
    // Preconditioning:
    //
    //     M z_k = r_k
    //
    // Search direction:
    //
    //     d_0 = z_0
    //
    // Pressure and residual updates:
    //
    //     alpha_k = (r_k,z_k) / (d_k,A d_k)
    //
    //     p_(k+1) = p_k + alpha_k d_k
    //
    //     r_(k+1) = r_k - alpha_k A d_k
    //
    // Direction update:
    //
    //     beta_k = (r_(k+1),z_(k+1)) / (r_k,z_k)
    //
    //     d_(k+1) = z_(k+1) + beta_k d_k
    //
    // Inputs:
    //     pres        initial pressure guess
    //     rhs1        pressure right-hand side
    //     pdiag       pressure diagonal and Jacobi preconditioner
    //     gstif       compact pressure off-diagonal coefficients
    //     bc_list     prescribed/reference pressure nodes
    //     bc_values   prescribed pressure values
    //
    // Output:
    //     pres        converged pressure field
    //
    // The returned Result contains Pressure CG iteration counts, residuals and
    // detailed timing information.
    //=========================================================================
    ConjugateGradient::Result ConjugateGradient::solvePressure(CBSStateSI& s)
    {
        const auto total_start = std::chrono::steady_clock::now();

        validate_cg_config(s);

        Result result;

        const auto setup_start = std::chrono::steady_clock::now();

        std::vector<char> active;
        std::vector<char> fixed;
        Array1D<Real> fixed_values;

        // Construct the active fluid-pressure space and the prescribed
        // pressure constraints.
        build_pressure_active_mask(s, active);
        build_fixed_pressure_mask(s, active, fixed, fixed_values);

        Array1D<Real> r;
        Array1D<Real> z;
        Array1D<Real> p;
        Array1D<Real> Ap;

        r.resize(s.cfg.npoin);
        z.resize(s.cfg.npoin);
        p.resize(s.cfg.npoin);
        Ap.resize(s.cfg.npoin);

        result.setup_seconds += elapsed_seconds(
            setup_start,
            std::chrono::steady_clock::now());

        {
            const auto t0 = std::chrono::steady_clock::now();
            // Apply solid-node inactivity and prescribed-pressure rows.
            impose_pressure_constraints(s, active, fixed, fixed_values);
            result.constraint_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
            // Calculate A p_0 for the initial residual.
            MatrixVectorCalc::pressureMultiply(s, s.pres, Ap);
            result.matvec_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            ++result.matvec_calls;
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                // Initial residual:
                //
                //     r_0 = b - A p_0
                r(ip) = is_pressure_active(active, ip) ? s.rhs1(ip) - Ap(ip) : 0.0;
            }

            zero_inactive_or_fixed_nodes(r, active, fixed, s.cfg.npoin);
            zero_inactive_or_fixed_nodes(Ap, active, fixed, s.cfg.npoin);
            result.vector_update_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
            result.initial_l2 = l2_norm_serial(r, active, fixed, s.cfg.npoin);
            result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            ++result.dot_calls;
        }

        if (result.initial_l2 <= 1.0e-300)
        {
            result.final_l2 = result.initial_l2;
            result.final_relative_l2 = 0.0;

            {
                const auto t0 = std::chrono::steady_clock::now();
                result.final_max_abs = max_abs_serial(r, active, fixed, s.cfg.npoin);
                result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                ++result.dot_calls;
            }

            result.converged = true;
            result.total_seconds = elapsed_seconds(total_start, std::chrono::steady_clock::now());
            return result;
        }

        result.final_l2 = result.initial_l2;
        result.final_relative_l2 = 1.0;

        {
            const auto t0 = std::chrono::steady_clock::now();
            result.final_max_abs = max_abs_serial(r, active, fixed, s.cfg.npoin);
            result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            ++result.dot_calls;
        }

        if (convergence_reached(s, result.final_relative_l2, result.final_max_abs))
        {
            result.converged = true;
            result.total_seconds = elapsed_seconds(total_start, std::chrono::steady_clock::now());
            return result;
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
            // Solve M z_0 = r_0.
            apply_preconditioner(s, r, z, active, fixed);
            result.preconditioner_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            ++result.preconditioner_calls;
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                // Initial CG search direction:
                //
                //     d_0 = z_0
                p(ip) = z(ip);
            }
            zero_inactive_or_fixed_nodes(p, active, fixed, s.cfg.npoin);
            result.vector_update_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
        }

        Real rz_old = 0.0;
        {
            const auto t0 = std::chrono::steady_clock::now();
            rz_old = dot_owned_serial(r, z, active, fixed, s.cfg.npoin);
            result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            ++result.dot_calls;
        }

        if (rz_old <= 1.0e-300 || !std::isfinite(rz_old))
        {
            throw std::runtime_error(
                "ConjugateGradient::solvePressure - invalid initial preconditioned residual dot product. "
                + diagnose_initial_cg_state(s, r, z, active, fixed, rz_old));
        }

        for (Int iter = 1; iter <= s.cfg.cg_max_iter; ++iter)
        {
            {
                const auto t0 = std::chrono::steady_clock::now();
                // Calculate A d_k.
                MatrixVectorCalc::pressureMultiply(s, p, Ap);
                result.matvec_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                ++result.matvec_calls;
            }

            {
                const auto t0 = std::chrono::steady_clock::now();
                zero_inactive_or_fixed_nodes(Ap, active, fixed, s.cfg.npoin);
                result.vector_update_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            }

            Real pAp = 0.0;
            {
                const auto t0 = std::chrono::steady_clock::now();
                pAp = dot_owned_serial(p, Ap, active, fixed, s.cfg.npoin);
                result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                ++result.dot_calls;
            }

            if (std::abs(pAp) <= 1.0e-300 || !std::isfinite(pAp))
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - breakdown: invalid pAp");
            }

            // CG step length:
            //
            //     alpha_k = (r_k,z_k) / (d_k,A d_k)
            const Real alpha = rz_old / pAp;

            {
                const auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (is_free_pressure_node(active, fixed, ip))
                    {
                        // Update pressure and residual:
                        //
                        //     p_(k+1) = p_k + alpha_k d_k
                        //
                        //     r_(k+1) = r_k - alpha_k A d_k
                        s.pres(ip) += alpha * p(ip);
                        r(ip) -= alpha * Ap(ip);
                    }
                }
                result.vector_update_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            }

            {
                const auto t0 = std::chrono::steady_clock::now();
                impose_pressure_constraints(s, active, fixed, fixed_values);
                zero_inactive_or_fixed_nodes(r, active, fixed, s.cfg.npoin);
                result.constraint_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            }

            result.iterations = iter;

            {
                const auto t0 = std::chrono::steady_clock::now();
                result.final_l2 = l2_norm_serial(r, active, fixed, s.cfg.npoin);
                result.final_max_abs = max_abs_serial(r, active, fixed, s.cfg.npoin);
                result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                result.dot_calls += 2;
            }

            result.final_relative_l2 = result.final_l2 / result.initial_l2;

            if (convergence_reached(s, result.final_relative_l2, result.final_max_abs))
            {
                result.converged = true;
                break;
            }

            {
                const auto t0 = std::chrono::steady_clock::now();
                apply_preconditioner(s, r, z, active, fixed);
                result.preconditioner_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                ++result.preconditioner_calls;
            }

            Real rz_new = 0.0;
            {
                const auto t0 = std::chrono::steady_clock::now();
                rz_new = dot_owned_serial(r, z, active, fixed, s.cfg.npoin);
                result.dot_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
                ++result.dot_calls;
            }

            if (rz_new <= 1.0e-300 || !std::isfinite(rz_new))
            {
                throw std::runtime_error(
                    "ConjugateGradient::solvePressure - invalid updated preconditioned residual dot product. "
                    + diagnose_initial_cg_state(s, r, z, active, fixed, rz_new));
            }

            // CG direction coefficient:
            //
            //     beta_k =
            //       (r_(k+1),z_(k+1)) / (r_k,z_k)
            const Real beta = rz_new / rz_old;

            {
                const auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (is_free_pressure_node(active, fixed, ip))
                    {
                        // New search direction:
                        //
                        //     d_(k+1) = z_(k+1) + beta_k d_k
                        p(ip) = z(ip) + beta * p(ip);
                    }
                    else
                    {
                        p(ip) = 0.0;
                    }
                }
                result.vector_update_seconds += elapsed_seconds(t0, std::chrono::steady_clock::now());
            }

            rz_old = rz_new;
        }

        result.total_seconds = elapsed_seconds(total_start, std::chrono::steady_clock::now());

        if (!result.converged)
        {
            std::cout << "WARNING: Pressure CG did not converge in "
                      << s.cfg.cg_max_iter
                      << " iterations. Final relative L2 = "
                      << result.final_relative_l2
                      << ", final max abs = "
                      << result.final_max_abs
                      << "\n";
        }

        return result;
    }
}
