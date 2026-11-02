// CBS3D++ G5a: isolated nonzero-velocity, multi-element thermal transport
// benchmark using the production EnergyAssembly::elementOperator and RHS.
// No CBS momentum/pressure solve and no production discretisation changes.
//
// Steady constant-coefficient convection-diffusion on [0,1]x[0,W]^2:
//   rhoCp * U * dT/dx = k * d^2T/dx^2
//   T(0)=1, T(1)=2, insulated transverse walls, U=(0.3,0,0).
// Exact: T(x) = 1 + expm1(Pe*x)/expm1(Pe), Pe=rhoCp*U/k=3.
//
// For each conformal tetra mesh, directly solve the assembled steady
// residual with an independent pivoted dense linear solver, and then verify:
//  1. actual production RHS vanishes on all unconstrained nodes;
//  2. uniform temperature is invariant for imposed divergence-free U;
//  3. solution obeys physical temperature bounds;
//  4. maximum nodal analytic error decreases under x-refinement;
//  5. nonzero CBS characteristic term (tau=dt/2), with dt proportional to h.
// This is a stationary prescribed-flow test, not an end-to-end coupled run.
// The added streamline diffusion is O(h); no order-2 claim is made.
#include "cbs/assembly/EnergyAssembly.hpp"
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
    using cbs::Int;
    using cbs::Real;

    constexpr Real u0 = 0.3;
    constexpr Real conductivity = 0.1;
    constexpr Real rhocp = 1.0;
    constexpr Real width = 0.25;
    constexpr Real temperature_left = 1.0;
    constexpr Real temperature_right = 2.0;
    constexpr Real courant = 0.3;

    struct Point
    {
        Real x, y, z;
    };

    using Tet = std::array<Int,4>;
    using Vec = std::array<Real,3>;

    Vec subtract(const Vec& a, const Vec& b)
    {
        return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
    }

    Vec cross(const Vec& a, const Vec& b)
    {
        return {{a[1]*b[2]-a[2]*b[1],
                 a[2]*b[0]-a[0]*b[2],
                 a[0]*b[1]-a[1]*b[0]}};
    }

    Real dot(const Vec& a, const Vec& b)
    {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    }

    Vec as_vector(const Point& p)
    {
        return {{p.x,p.y,p.z}};
    }

    Real exact_temperature(const Real x)
    {
        const Real pe = rhocp*u0/conductivity;
        return temperature_left+
            (temperature_right-temperature_left)*
            std::expm1(pe*x)/std::expm1(pe);
    }

    struct Mesh
    {
        std::vector<Point> points;
        std::vector<Tet> elements;
    };

    Mesh make_mesh(const Int nx)
    {
        constexpr Int ny=2, nz=2;
        Mesh mesh;
        const auto id = [&](const Int i, const Int j, const Int k) -> Int {
            return 1 + (i*(ny+1)+j)*(nz+1)+k;
        };
        for (Int i=0;i<=nx;++i)
        {
            for (Int j=0;j<=ny;++j)
            {
                for (Int k=0;k<=nz;++k)
                {
                    mesh.points.push_back({
                        Real(i)/Real(nx),
                        width*Real(j)/Real(ny),
                        width*Real(k)/Real(nz)
                    });
                }
            }
        }
        for(Int i=0;i<nx;++i)
        {
            for(Int j=0;j<ny;++j)
            {
                for(Int k=0;k<nz;++k)
                {
                    const Int a=id(i,j,k),b=id(i+1,j,k),
                        c=id(i+1,j+1,k),d=id(i,j+1,k),
                        e=id(i,j,k+1),f=id(i+1,j,k+1),
                        g=id(i+1,j+1,k+1),h=id(i,j+1,k+1);
                    const std::array<Tet,6> tets = {{
                        {{a,b,c,g}},{{a,c,d,g}},{{a,d,h,g}},
                        {{a,h,e,g}},{{a,e,f,g}},{{a,f,b,g}}
                    }};
                    mesh.elements.insert(mesh.elements.end(),
                        tets.begin(),tets.end());
                }
            }
        }
        return mesh;
    }

    CBSStateSI prepare_state(const Mesh& mesh, const Int nx)
    {
        CBSStateSI s;
        s.set_problem_sizes(
            static_cast<Int>(mesh.elements.size()),
            static_cast<Int>(mesh.points.size()),0,0);
        s.cfg.temp_calc=1;
        s.cfg.source_solid=0.0;
        s.cfg.heat_flux_bc=0.0;
        s.cfg.turbulence_on=0;
        s.cfg.turbulent_thermal_diffusivity_on=0;
        s.cfg.inlet_temperature=temperature_left;
        s.unkno.fill(0.0);
        s.temperature.fill(temperature_left);
        s.temperature1.fill(temperature_left);
        s.rhs1.fill(0.0);
        for(Int ip=1;ip<=s.cfg.npoin;++ip)
        {
            s.unkno(1,ip)=u0;
        }

        for(Int ie=1;ie<=s.cfg.nelem;++ie)
        {
            const Tet& tet=mesh.elements[static_cast<std::size_t>(ie-1)];
            const Vec p0=as_vector(mesh.points[static_cast<std::size_t>(tet[0]-1)]);
            const Vec a=subtract(as_vector(mesh.points[static_cast<std::size_t>(tet[1]-1)]),p0);
            const Vec b=subtract(as_vector(mesh.points[static_cast<std::size_t>(tet[2]-1)]),p0);
            const Vec c=subtract(as_vector(mesh.points[static_cast<std::size_t>(tet[3]-1)]),p0);
            const Real det=dot(a,cross(b,c));
            if(std::abs(det)<1.0e-16)
            {
                throw std::runtime_error("G5a degenerate tetrahedron");
            }
            const Vec g1={{cross(b,c)[0]/det,cross(b,c)[1]/det,cross(b,c)[2]/det}};
            const Vec g2={{cross(c,a)[0]/det,cross(c,a)[1]/det,cross(c,a)[2]/det}};
            const Vec g3={{cross(a,b)[0]/det,cross(a,b)[1]/det,cross(a,b)[2]/det}};
            const Vec g0={{-g1[0]-g2[0]-g3[0],
                           -g1[1]-g2[1]-g3[1],
                           -g1[2]-g2[2]-g3[2]}};
            const std::array<Vec,4> grads={{g0,g1,g2,g3}};
            for(Int local=1;local<=4;++local)
            {
                s.intma(local,ie)=tet[static_cast<std::size_t>(local-1)];
                for(Int dimension=1;dimension<=3;++dimension)
                {
                    s.dNkdx((ie-1)*12+(dimension-1)*4+local)=
                        grads[static_cast<std::size_t>(local-1)]
                            [static_cast<std::size_t>(dimension-1)];
                }
            }
            s.detJ(ie)=std::abs(det);
            s.delte(ie)=courant/(Real(nx)*u0);
            s.mat_elem(ie)=0;
            s.rho_e(ie)=1.0;
            s.cp_e(ie)=rhocp;
            s.rho_cp_e(ie)=rhocp;
            s.k_e(ie)=conductivity;
            s.k_eff_e(ie)=conductivity;
            s.Qvol_e(ie)=0.0;
        }
        return s;
    }

    // Test-only dense Gaussian elimination; production solver linear algebra
    // is deliberately excluded.  Each row is A=-dR/dT; Dirichlet x-faces
    // are applied by row replacement.
    std::vector<Real> solve_system(
        std::vector<Real> matrix,std::vector<Real> rhs,const Int n)
    {
        for(Int i=0;i<n;++i)
        {
            Int pivot=i;
            for(Int j=i+1;j<n;++j)
            {
                if(std::abs(matrix[static_cast<std::size_t>(j)*n+i])>
                   std::abs(matrix[static_cast<std::size_t>(pivot)*n+i]))
                {
                    pivot=j;
                }
            }
            if(std::abs(matrix[static_cast<std::size_t>(pivot)*n+i])<1.0e-14)
            {
                throw std::runtime_error("G5a assembled matrix singular");
            }
            if(pivot!=i)
            {
                for(Int j=i;j<n;++j)
                    std::swap(matrix[static_cast<std::size_t>(i)*n+j],
                              matrix[static_cast<std::size_t>(pivot)*n+j]);
                std::swap(rhs[static_cast<std::size_t>(i)],rhs[static_cast<std::size_t>(pivot)]);
            }
            for(Int j=i+1;j<n;++j)
            {
                const Real factor=
                    matrix[static_cast<std::size_t>(j)*n+i]/
                    matrix[static_cast<std::size_t>(i)*n+i];
                if(factor==0.0)continue;
                matrix[static_cast<std::size_t>(j)*n+i]=0.0;
                for(Int k=i+1;k<n;++k)
                {
                    matrix[static_cast<std::size_t>(j)*n+k]-=
                        factor*matrix[static_cast<std::size_t>(i)*n+k];
                }
                rhs[static_cast<std::size_t>(j)]-=
                    factor*rhs[static_cast<std::size_t>(i)];
            }
        }
        std::vector<Real> x(static_cast<std::size_t>(n),0.0);
        for(Int i=n-1;i>=0;--i)
        {
            Real remaining=rhs[static_cast<std::size_t>(i)];
            for(Int j=i+1;j<n;++j)
            {
                remaining-=matrix[static_cast<std::size_t>(i)*n+j]*
                    x[static_cast<std::size_t>(j)];
            }
            x[static_cast<std::size_t>(i)]=
                remaining/matrix[static_cast<std::size_t>(i)*n+i];
        }
        return x;
    }

    struct Verification
    {
        Real max_error=0.0;
        Real max_residual=0.0;
        Real constant_residual=0.0;
        Real min_temperature=std::numeric_limits<Real>::infinity();
        Real max_temperature=-std::numeric_limits<Real>::infinity();
        Int npoints=0;
        Int ntets=0;
    };

    Verification run(const Int nx)
    {
        const Mesh mesh=make_mesh(nx);
        CBSStateSI s=prepare_state(mesh,nx);
        const Int n=s.cfg.npoin;
        std::vector<Real> A(static_cast<std::size_t>(n)*n,0.0);
        std::vector<Real> rhs(static_cast<std::size_t>(n),0.0);
        std::vector<char> fixed(static_cast<std::size_t>(n),0);

        for(Int ie=1;ie<=s.cfg.nelem;++ie)
        {
            Real M[5][5]={{0.0}};
            cbs::EnergyAssembly::elementOperator(s,ie,M);
            for(Int a=1;a<=4;++a)
            {
                const Int i=s.intma(a,ie)-1;
                for(Int b=1;b<=4;++b)
                {
                    const Int j=s.intma(b,ie)-1;
                    A[static_cast<std::size_t>(i)*n+j]-=M[a][b];
                }
            }
        }

        // Independent constant-state check with nonzero u.
        s.temperature1.fill(temperature_left);
        cbs::EnergyAssembly::assembleStep4Rhs(s);
        Verification result;
        result.npoints=n;
        result.ntets=s.cfg.nelem;
        for(Int ip=1;ip<=n;++ip)
        {
            result.constant_residual=std::max(result.constant_residual,
                                             std::abs(s.rhs1(ip)));
        }

        for(Int ip=0;ip<n;++ip)
        {
            const Real x=mesh.points[static_cast<std::size_t>(ip)].x;
            if(x==0.0 || x==1.0)
            {
                fixed[static_cast<std::size_t>(ip)]=1;
                for(Int j=0;j<n;++j)
                    A[static_cast<std::size_t>(ip)*n+j]=0.0;
                A[static_cast<std::size_t>(ip)*n+ip]=1.0;
                rhs[static_cast<std::size_t>(ip)]=
                    x==0.0?temperature_left:temperature_right;
            }
        }
        const auto solution=solve_system(A,rhs,n);
        for(Int ip=1;ip<=n;++ip)
        {
            const Real t=solution[static_cast<std::size_t>(ip-1)];
            const Real x=mesh.points[static_cast<std::size_t>(ip-1)].x;
            s.temperature1(ip)=t;
            result.max_error=std::max(result.max_error,std::abs(t-exact_temperature(x)));
            result.min_temperature=std::min(result.min_temperature,t);
            result.max_temperature=std::max(result.max_temperature,t);
        }

        // Verify the production matrix agrees with its independently
        // executed residual on all unconstrained nodes.
        cbs::EnergyAssembly::assembleStep4Rhs(s);
        for(Int ip=1;ip<=n;++ip)
        {
            if(!fixed[static_cast<std::size_t>(ip-1)])
            {
                result.max_residual=std::max(result.max_residual,
                                             std::abs(s.rhs1(ip)));
            }
        }
        std::cout<<std::setprecision(11)<<std::scientific
                 <<"nx="<<nx<<" nodes="<<n<<" tetrahedra="<<s.cfg.nelem
                 <<" max_error_K="<<result.max_error
                 <<" residual_max="<<result.max_residual
                 <<" uniform_residual_max="<<result.constant_residual
                 <<" Tmin="<<result.min_temperature
                 <<" Tmax="<<result.max_temperature<<"\n";
        return result;
    }
}

int main()
{
    try
    {
        const Verification coarse=run(8);
        const Verification medium=run(16);
        const Verification fine=run(32);
        const bool converged=
            coarse.max_error>1.3*medium.max_error &&
            medium.max_error>1.4*fine.max_error &&
            fine.max_error<7.5e-3;
        const bool residual_ok=
            std::max({coarse.max_residual,medium.max_residual,fine.max_residual})<1.0e-10 &&
            std::max({coarse.constant_residual,medium.constant_residual,
                      fine.constant_residual})<1.0e-11;
        const bool bounded=
            std::min({coarse.min_temperature,medium.min_temperature,fine.min_temperature})>=temperature_left-1.0e-10 &&
            std::max({coarse.max_temperature,medium.max_temperature,fine.max_temperature})<=temperature_right+1.0e-10;
        std::cout<<"Mesh-error ratio 8/16 = "<<coarse.max_error/medium.max_error
                 <<" 16/32 = "<<medium.max_error/fine.max_error<<"\n";
        std::cout<<"Analytic spatial refinement: "<<(converged?"PASS":"FAIL")<<"\n";
        std::cout<<"Production residual consistency: "<<(residual_ok?"PASS":"FAIL")<<"\n";
        std::cout<<"Convection-diffusion bounds (Pe=3): "<<(bounded?"PASS":"FAIL")<<"\n";
        if(!(converged && residual_ok && bounded))
        {
            std::cerr<<"G5A PRESCRIBED-VELOCITY ADVECTION-DIFFUSION FAIL\n";
            return 1;
        }
        std::cout<<"G5A PRESCRIBED-VELOCITY ADVECTION-DIFFUSION PASS\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"G5a ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
