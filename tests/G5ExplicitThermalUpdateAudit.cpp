// CBS3D++ verification G5b: explicit first-step positivity audit.
//
// This isolated test calls the UNMODIFIED production Step-4 thermal
// elementOperator / assembleStep4Rhs and TimeStep::updateLhsDiagonal.
// A test-only P1 tetrahedral fluid mesh has nonzero prescribed velocity.
// The experiment identifies negative off-diagonal coefficients in the
// full first-step matrix B = I + diag(elcoe2p) * dR/dT.
//
// If an unconstrained B(i,j)<0, an admissible temperature hat function
// T_j=T_ref+25K and T_{k!=j}=T_ref gives T_i^{new}<T_ref.
// Such a result shows the discrete maximum principle is not guaranteed
// by the *unlimited* explicit operator for that mesh/Courant/Pe setup.
// It does NOT imply unbounded growth, instability, or that Job 29706
// has this mechanism; those require separate evidence.
//
// PASS means matrix/direct-update consistency and correct capacities.
// It deliberately does NOT mean the method is positivity preserving.

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/timestep/TimeStep.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Real;
    using cbs::Int;
    using Point = std::array<Real,3>;
    using Tet = std::array<Int,4>;

    constexpr Real reference_temperature=298.15;
    constexpr Real bump_temperature=25.0;
    constexpr Real speed=0.3;
    constexpr Real rhocp=1.0;
    constexpr Real width=0.25;

    Point minus(const Point& a,const Point& b)
    {
        return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
    }
    Point cross(const Point& a,const Point& b)
    {
        return {{a[1]*b[2]-a[2]*b[1],
                 a[2]*b[0]-a[0]*b[2],
                 a[0]*b[1]-a[1]*b[0]}};
    }
    Real dot(const Point& a,const Point& b)
    {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    }

    struct Mesh
    {
        std::vector<Point> points;
        std::vector<Tet> cells;
    };

    Mesh mesh(const Int nx)
    {
        constexpr Int ny=2,nz=2;
        Mesh result;
        const auto id=[](const Int i,const Int j,const Int k)->Int
        {
            return 1+(i*(ny+1)+j)*(nz+1)+k;
        };

        for(Int i=0;i<=nx;++i)
            for(Int j=0;j<=ny;++j)
                for(Int k=0;k<=nz;++k)
                    result.points.push_back({{
                        Real(i)/nx, width*Real(j)/ny,width*Real(k)/nz
                    }});

        for(Int i=0;i<nx;++i)
            for(Int j=0;j<ny;++j)
                for(Int k=0;k<nz;++k)
                {
                    const Int a=id(i,j,k),b=id(i+1,j,k),
                              c=id(i+1,j+1,k),d=id(i,j+1,k),
                              e=id(i,j,k+1),f=id(i+1,j,k+1),
                              g=id(i+1,j+1,k+1),h=id(i,j+1,k+1);
                    const std::array<Tet,6> t={{
                        {{a,b,c,g}},{{a,c,d,g}},{{a,d,h,g}},
                        {{a,h,e,g}},{{a,e,f,g}},{{a,f,b,g}}
                    }};
                    result.cells.insert(result.cells.end(),t.begin(),t.end());
                }
        return result;
    }

    bool interior_x(const Point& p)
    {
        return p[0]>1.0e-12 && p[0]<1.0-1.0e-12;
    }

    CBSStateSI state(const Mesh& m,const Int nx,
                     const Real k,const Real courant,
                     std::vector<Real>& reference_capacity)
    {
        CBSStateSI s;
        s.set_problem_sizes(
            static_cast<Int>(m.cells.size()),
            static_cast<Int>(m.points.size()),0,0);
        s.cfg.temp_calc=1;
        s.cfg.bc_fixed=0;
        s.cfg.source_solid=0.0;
        s.cfg.heat_flux_bc=0.0;
        s.cfg.turbulence_on=0;
        s.cfg.turbulent_thermal_diffusivity_on=0;
        s.cfg.inlet_temperature=reference_temperature;

        s.unkno.fill(0.0);
        s.temperature.fill(reference_temperature);
        s.temperature1.fill(reference_temperature);
        s.pres.fill(0.0);
        s.rhs1.fill(0.0);
        s.gstifE.fill(0.0);
        s.pdiagE.fill(0.0);

        for(Int ip=1;ip<=s.cfg.npoin;++ip)
            s.unkno(1,ip)=speed;

        const Real dt=courant/(nx*speed);
        reference_capacity.assign(static_cast<std::size_t>(s.cfg.npoin),0.0);

        for(Int ie=1;ie<=s.cfg.nelem;++ie)
        {
            const Tet& tet=m.cells[static_cast<std::size_t>(ie-1)];
            const Point p0=m.points[static_cast<std::size_t>(tet[0]-1)];
            const Point a=minus(m.points[static_cast<std::size_t>(tet[1]-1)],p0);
            const Point b=minus(m.points[static_cast<std::size_t>(tet[2]-1)],p0);
            const Point c=minus(m.points[static_cast<std::size_t>(tet[3]-1)],p0);
            const Real d=dot(a,cross(b,c));
            if(std::abs(d)<1.0e-16)
                throw std::runtime_error("degenerate test tetrahedron");

            const Point ca=cross(b,c),cb=cross(c,a),cc=cross(a,b);
            const Point g1={{ca[0]/d,ca[1]/d,ca[2]/d}};
            const Point g2={{cb[0]/d,cb[1]/d,cb[2]/d}};
            const Point g3={{cc[0]/d,cc[1]/d,cc[2]/d}};
            const Point g0={{-g1[0]-g2[0]-g3[0],
                             -g1[1]-g2[1]-g3[1],
                             -g1[2]-g2[2]-g3[2]}};
            const std::array<Point,4> gradients={{g0,g1,g2,g3}};

            s.detJ(ie)=std::abs(d);
            s.delte(ie)=dt;
            s.mat_elem(ie)=0;
            s.rho_e(ie)=1.0;
            s.cp_e(ie)=rhocp;
            s.rho_cp_e(ie)=rhocp;
            s.k_e(ie)=k;
            s.k_eff_e(ie)=k;
            s.Qvol_e(ie)=0.0;

            for(Int a_local=1;a_local<=4;++a_local)
            {
                const Int ip=tet[static_cast<std::size_t>(a_local-1)];
                s.intma(a_local,ie)=ip;
                for(Int dim=1;dim<=3;++dim)
                    s.dNkdx((ie-1)*12+(dim-1)*4+a_local)=
                        gradients[static_cast<std::size_t>(a_local-1)]
                                 [static_cast<std::size_t>(dim-1)];

                const Real lumped=std::abs(d)/24.0;
                s.elcoe_e((ie-1)*4+a_local)=lumped;
                // Pressure diagonals are nonzero only to make the production
                // mass/pressure bookkeeping well-defined for fluid nodes.
                s.pdiagE((ie-1)*4+a_local)=1.0;
                reference_capacity[static_cast<std::size_t>(ip-1)]+=
                    rhocp*lumped/dt;
            }
        }
        return s;
    }

    bool run_case(const Int nx,const Real k,const Real courant)
    {
        const Mesh m=mesh(nx);
        std::vector<Real> capacity;
        CBSStateSI s=state(m,nx,k,courant,capacity);
        cbs::TimeStep::updateLhsDiagonal(s);

        bool capacity_ok=true;
        Real max_capacity_error=0.0;
        for(Int ip=1;ip<=s.cfg.npoin;++ip)
        {
            const Real reference_inverse=1.0/capacity[static_cast<std::size_t>(ip-1)];
            max_capacity_error=std::max(
                max_capacity_error,std::abs(s.elcoe2p(ip)-reference_inverse));
            capacity_ok &= std::abs(s.elcoe2p(ip)-reference_inverse)<
                1.0e-11*std::max(1.0,std::abs(reference_inverse));
        }

        // Collect the actual production residual Jacobian, with the sign
        // convention dR/dT (negative stiffness).
        const Int n=s.cfg.npoin;
        std::vector<Real> matrix(static_cast<std::size_t>(n)*n,0.0);
        for(Int ie=1;ie<=s.cfg.nelem;++ie)
        {
            Real element[5][5]={{0.0}};
            cbs::EnergyAssembly::elementOperator(s,ie,element);
            for(Int a=1;a<=4;++a)
            {
                const Int ip=s.intma(a,ie)-1;
                for(Int b=1;b<=4;++b)
                {
                    const Int jp=s.intma(b,ie)-1;
                    matrix[static_cast<std::size_t>(ip)*n+jp]+=element[a][b];
                }
            }
        }

        s.temperature1.fill(reference_temperature);
        cbs::EnergyAssembly::assembleStep4Rhs(s);
        Real uniform_rhs=0.0;
        for(Int ip=1;ip<=n;++ip)
            uniform_rhs=std::max(uniform_rhs,std::abs(s.rhs1(ip)));

        // Monotonicity of the FULL explicit update on interior degrees
        // requires all its off-diagonal entries to be nonnegative.
        // Report negative entries, select the strongest as a counterexample.
        Int negative_count=0,best_i=-1,best_j=-1;
        Real most_negative=0.0;
        for(Int ip=0;ip<n;++ip)
        {
            if(!interior_x(m.points[static_cast<std::size_t>(ip)]))
                continue;
            for(Int jp=0;jp<n;++jp)
            {
                if(ip==jp ||
                   !interior_x(m.points[static_cast<std::size_t>(jp)]))
                    continue;
                const Real off=s.elcoe2p(ip+1)*
                    matrix[static_cast<std::size_t>(ip)*n+jp];
                if(off < -1.0e-13)
                {
                    ++negative_count;
                    if(off < most_negative)
                    {
                        most_negative=off;
                        best_i=ip;
                        best_j=jp;
                    }
                }
            }
        }

        Real observed=reference_temperature;
        Real predicted=reference_temperature;
        Real residual_prediction_error=0.0;
        if(best_i>=0)
        {
            // Bounded initial P1 hat function: one hot interior node, all
            // others equal to Tref, no negative heat sources or external flux.
            s.temperature1.fill(reference_temperature);
            s.temperature1(best_j+1)=reference_temperature+bump_temperature;
            cbs::EnergyAssembly::assembleStep4Rhs(s);
            observed=s.temperature1(best_i+1)+
                s.elcoe2p(best_i+1)*s.rhs1(best_i+1);
            predicted=reference_temperature+bump_temperature*most_negative;
            residual_prediction_error=std::abs(observed-predicted);
        }

        const Real pe=rhocp*speed/k;
        std::cout<<std::scientific<<std::setprecision(11)
                 <<"nx="<<nx<<" Pe_L="<<pe
                 <<" Courant_x="<<courant
                 <<" dt="<<s.delte(1)
                 <<" max_cap_error="<<max_capacity_error
                 <<" uniform_rhs="<<uniform_rhs
                 <<" negative_offdiag_count="<<negative_count;
        if(best_i>=0)
        {
            std::cout<<" most_negative_Bij="<<most_negative
                     <<" i="<<best_i+1<<" j="<<best_j+1
                     <<" T_i_after_hat="<<observed
                     <<" T_i_predicted="<<predicted
                     <<" analytic_linearization_error="<<residual_prediction_error
                     <<" Tmin_bound_violation_K="<<
                        std::max(0.0,reference_temperature-observed);
        }
        else
        {
            std::cout<<" no_counterexample_for_selected_parameters";
        }
        std::cout<<"\n";

        const bool ok=capacity_ok && uniform_rhs<1.0e-10 &&
            residual_prediction_error<1.0e-9 &&
            (best_i<0 || observed<reference_temperature-1.0e-12);
        return ok;
    }
}

int main()
{
    try
    {
        bool ok=true;
        for(const Int nx : {8,16})
            for(const Real k : {0.1,0.001})
                for(const Real courant : {0.1,0.3,0.6})
                    ok=run_case(nx,k,courant)&&ok;
        if(!ok)
        {
            std::cerr<<"G5B PRODUCTION EXPLICIT UPDATE AUDIT FAIL\n";
            return 1;
        }
        std::cout<<"G5B PRODUCTION EXPLICIT UPDATE AUDIT PASS\n"
                 <<"PASS does NOT assert DMP: inspect negative_offdiag_count "
                    "and first-step hat-perturbation undershoot above.\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"G5b ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
