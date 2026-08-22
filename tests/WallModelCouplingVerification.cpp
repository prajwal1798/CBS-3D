//=============================================================================
// Production wall-coupling regression.
//
// Verifies three independent contracts:
//   1. a modelled TRI3 replaces, rather than adds to, the existing natural
//      viscous face contribution;
//   2. the three P1 nodal loads sum exactly to the modelled face load and the
//      traction is dissipative;
//   3. impermeability is the projector onto the span of incident wall normals,
//      leaving tangent-space velocity free on planes/edges and zero only at a
//      rank-three corner.
//=============================================================================

#include "cbs/turbulence/WallModelCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;
    using cbs::turbulence::WallModelCoupling;

    bool close(
        const Real a,
        const Real b,
        const Real relative = 2.0e-12,
        const Real absolute = 2.0e-13)
    {
        return std::abs(a - b) <=
            absolute + relative * std::max(std::abs(a), std::abs(b));
    }

    void set_standard_tet_face_map(CBSStateSI& s)
    {
        s.ippn1.resize(4, 3);

        s.ippn1(1, 1) = 2;
        s.ippn1(1, 2) = 3;
        s.ippn1(1, 3) = 4;

        s.ippn1(2, 1) = 1;
        s.ippn1(2, 2) = 4;
        s.ippn1(2, 3) = 3;

        s.ippn1(3, 1) = 1;
        s.ippn1(3, 2) = 2;
        s.ippn1(3, 3) = 4;

        s.ippn1(4, 1) = 1;
        s.ippn1(4, 2) = 3;
        s.ippn1(4, 3) = 2;
    }

    Int gradient_index(
        const CBSStateSI& s,
        const Int dim,
        const Int local_node)
    {
        return (dim - 1) * s.cfg.nep + local_node;
    }

    CBSStateSI base_state(const Int boundary_faces)
    {
        CBSStateSI s;

        s.cfg.npoin = 4;
        s.cfg.nelem = 1;
        s.cfg.nboun = boundary_faces;
        s.cfg.ndim = 3;
        s.cfg.nep = 4;
        s.cfg.nsid = 4;
        s.cfg.nsidp = 3;
        s.cfg.bsid = 6;
        s.cfg.nsidpl = 4;
        s.cfg.nsidpe = 5;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.dimensional_mode = 0;
        s.cfg.material_properties_enabled = 0;
        s.cfg.ani = 1.0e-5;

        s.intma.resize(4, 1);
        s.intma(1, 1) = 1;
        s.intma(2, 1) = 2;
        s.intma(3, 1) = 3;
        s.intma(4, 1) = 4;

        s.mat_elem.resize(1);
        s.mat_elem(1) = 0;

        s.coord.resize(3, 4);
        const std::array<std::array<Real, 3>, 4> xyz =
        {{
            {{0.0, 0.0, 0.0}},
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 0.2}}
        }};

        for (Int ip = 1; ip <= 4; ++ip)
        {
            for (Int dim = 1; dim <= 3; ++dim)
            {
                s.coord(dim, ip) =
                    xyz[static_cast<std::size_t>(ip - 1)]
                       [static_cast<std::size_t>(dim - 1)];
            }
        }

        set_standard_tet_face_map(s);

        s.detJ.resize(1);
        s.detJ(1) = 0.2;

        s.iside.resize(6, boundary_faces);
        s.face_norm.resize(4, boundary_faces);
        s.annxf.resize(4, 4, 1);
        s.annxf.fill(0.0);
        s.fedge.resize(4, 1);
        s.fedge.fill(0);

        s.dNkdx.resize(12);

        // Tetrahedron x=N2, y=N3, z=0.2 N4.
        const Real gradients[3][4] =
        {
            {-1.0, 1.0, 0.0, 0.0},
            {-1.0, 0.0, 1.0, 0.0},
            {-5.0, 0.0, 0.0, 5.0}
        };

        for (Int dim = 1; dim <= 3; ++dim)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                s.dNkdx(gradient_index(s, dim, a)) =
                    gradients[dim - 1][a - 1];
            }
        }

        s.unkn1.resize(3, 4);
        s.unkno.resize(3, 4);
        s.rhs.resize(3, 4);
        s.unkn1.fill(0.0);
        s.unkno.fill(0.0);
        s.rhs.fill(0.0);

        s.node_velocity_bc_type.resize(4);
        s.node_velocity_bc_type.fill(CBSStateSI::velocity_bc_free);

        return s;
    }

    void add_face(
        CBSStateSI& s,
        const Int ib,
        const Int local_face,
        const std::array<Int, 3>& nodes,
        const std::array<Real, 3>& unit_normal,
        const Real area)
    {
        for (Int i = 1; i <= 3; ++i)
        {
            s.iside(i, ib) = nodes[static_cast<std::size_t>(i - 1)];
            s.face_norm(i, ib) =
                area * unit_normal[static_cast<std::size_t>(i - 1)];
            s.annxf(i, local_face, 1) =
                area * unit_normal[static_cast<std::size_t>(i - 1)];
        }

        s.face_norm(4, ib) = area;
        s.annxf(4, local_face, 1) = area;
        s.iside(4, ib) = local_face;
        s.iside(5, ib) = 1;
        s.iside(6, ib) = s.cfg.bc_noslip_adiabatic_wall;
        s.fedge(local_face, 1) = 1;

        for (const Int ip : nodes)
        {
            s.node_velocity_bc_type(ip) = CBSStateSI::velocity_bc_noslip;
        }
    }

    bool verify_flux_replacement()
    {
        CBSStateSI s = base_state(1);

        // z=0 face, opposite local node 4.  The fluid occupies z>0, so the
        // outward normal is -ez and A=1/2.
        add_face(
            s,
            1,
            4,
            {1, 3, 2},
            {0.0, 0.0, -1.0},
            0.5);

        // u=10 at the first off-wall sample node and u=0 on the legacy wall.
        // Hence du/dz=50 and the OLD natural diffusion contribution at each
        // wall node is
        //
        //   nu (du/dz) (A n_z) / 3 = -8.333333...e-5.
        s.unkn1(1, 4) = 10.0;

        const Real old_natural =
            s.cfg.ani * 50.0 * (-0.5) / 3.0;

        for (const Int ip : {1, 2, 3})
        {
            s.rhs(1, ip) = old_natural;
        }

        const auto diagnostic =
            WallModelCoupling::replaceMomentumWallFlux(s);

        if (diagnostic.local_faces != 1 ||
            WallModelCoupling::globalWallFaceCount(s) != 1 ||
            !(diagnostic.modeled_wall_work <= 0.0) ||
            !(diagnostic.minimum_y_plus > 0.0) ||
            !(diagnostic.maximum_y_plus >= diagnostic.minimum_y_plus))
        {
            std::printf("FAIL: wall-face/model diagnostics\n");
            return false;
        }

        // If the old natural term had merely been left in place and the wall
        // load added on top, this equality would fail by exactly 3*old_natural.
        Real nodal_sum[3] = {0.0, 0.0, 0.0};

        for (const Int ip : {1, 2, 3})
        {
            for (Int dim = 1; dim <= 3; ++dim)
            {
                nodal_sum[dim - 1] += s.rhs(dim, ip);
            }
        }

        for (Int dim = 0; dim < 3; ++dim)
        {
            if (!close(
                    nodal_sum[dim],
                    diagnostic.modeled_surface_load[
                        static_cast<std::size_t>(dim)]))
            {
                std::printf(
                    "FAIL: face load conservation dim=%d sum=% .17e face=% .17e\n",
                    dim,
                    nodal_sum[dim],
                    diagnostic.modeled_surface_load[
                        static_cast<std::size_t>(dim)]);
                return false;
            }
        }

        for (const Int ip : {1, 2, 3})
        {
            if (!close(
                    s.rhs(1, ip),
                    diagnostic.modeled_surface_load[0] / 3.0) ||
                !close(s.rhs(2, ip), 0.0) ||
                !close(s.rhs(3, ip), 0.0))
            {
                std::printf("FAIL: constant TRI3 nodal traction\n");
                return false;
            }
        }

        if (!(diagnostic.modeled_surface_load[0] < 0.0))
        {
            std::printf("FAIL: wall traction does not oppose +x sample flow\n");
            return false;
        }

        return true;
    }

    bool verify_normal_span_projection(const bool flip_normals)
    {
        CBSStateSI s = base_state(3);
        const Real sign = flip_normals ? -1.0 : 1.0;

        // Three coordinate-plane walls meeting at node 1.
        add_face(s, 1, 4, {1, 3, 2}, {0.0, 0.0, -sign}, 0.5); // z=0
        add_face(s, 2, 3, {1, 2, 4}, {0.0, -sign, 0.0}, 0.1); // y=0
        add_face(s, 3, 2, {1, 4, 3}, {-sign, 0.0, 0.0}, 0.1); // x=0

        for (Int ip = 1; ip <= 4; ++ip)
        {
            s.unkno(1, ip) = 3.0;
            s.unkno(2, ip) = 4.0;
            s.unkno(3, ip) = 5.0;
        }

        const auto captured =
            WallModelCoupling::captureVelocity(s, false);

        // Simulate the existing strong no-slip package that the production
        // wrapper deliberately lets execute before restoring tangent velocity.
        for (Int ip = 1; ip <= 4; ++ip)
        {
            s.unkno(1, ip) = 0.0;
            s.unkno(2, ip) = 0.0;
            s.unkno(3, ip) = 0.0;
        }

        WallModelCoupling::restoreTangentialAndEnforceImpermeability(
            s,
            captured);

        // Node 1 sees three independent normals -> no tangent-space DOF.
        if (!close(s.unkno(1, 1), 0.0) ||
            !close(s.unkno(2, 1), 0.0) ||
            !close(s.unkno(3, 1), 0.0))
        {
            std::printf("FAIL: rank-three wall corner projection\n");
            return false;
        }

        // Node 2 lies on y=0 and z=0, so only x is free.
        if (!close(s.unkno(1, 2), 3.0) ||
            !close(s.unkno(2, 2), 0.0) ||
            !close(s.unkno(3, 2), 0.0))
        {
            std::printf("FAIL: two-normal wall-edge projection\n");
            return false;
        }

        // Node 3 lies on x=0 and z=0, so only y is free.
        if (!close(s.unkno(1, 3), 0.0) ||
            !close(s.unkno(2, 3), 4.0) ||
            !close(s.unkno(3, 3), 0.0))
        {
            std::printf("FAIL: second two-normal wall-edge projection\n");
            return false;
        }

        // Node 4 lies on x=0 and y=0, so only z is free.
        if (!close(s.unkno(1, 4), 0.0) ||
            !close(s.unkno(2, 4), 0.0) ||
            !close(s.unkno(3, 4), 5.0))
        {
            std::printf("FAIL: third two-normal wall-edge projection\n");
            return false;
        }

        return true;
    }

    bool verify_cht_rejected()
    {
        CBSStateSI s = base_state(1);
        add_face(s, 1, 4, {1, 3, 2}, {0.0, 0.0, -1.0}, 0.5);

        // The current production wall model is momentum-only and must not be
        // applied to CHT before thermal/interface treatment is validated.
        s.mat_elem(1) = 1;

        try
        {
            (void)WallModelCoupling::globalWallFaceCount(s);
        }
        catch (const std::exception&)
        {
            return true;
        }

        std::printf("FAIL: CHT wall-model activation was not rejected\n");
        return false;
    }
}

int main()
{
#if defined(_WIN32)
    _putenv_s("CBS3D_SA_WALL_TREATMENT", "1");
#else
    setenv("CBS3D_SA_WALL_TREATMENT", "1", 1);
#endif

    try
    {
        if (!verify_flux_replacement())
        {
            return 1;
        }

        if (!verify_normal_span_projection(false) ||
            !verify_normal_span_projection(true))
        {
            return 1;
        }

        if (!verify_cht_rejected())
        {
            return 1;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("FAIL: unexpected exception: %s\n", error.what());
        return 1;
    }

    std::printf("PASS: production SA wall-model coupling\n");
    return 0;
}
