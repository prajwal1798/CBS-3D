#pragma once

//=============================================================================
// CBS3D_SI
// Greedy element coloring for race-free OpenMP scatter assembly.
//
// Two elements that share any node must get different colors. Within a single
// color, every node is touched by at most one element, so the per-element
// scatter (rhs(.,ip) += ..., y(ip) += ...) can run in parallel with no data
// race and no atomics. Colors are processed sequentially; elements within a
// color run in parallel.
//
// Fills s.ncolor, s.color_ptr (size ncolor+1), s.color_elem (size nelem).
// Also fills pressure-only colour arrays containing only mat_elem == 0 elements.
// Build once after connectivity and material IDs are known.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <vector>

namespace cbs
{
    class Coloring
    {
    public:
        static void build(CBSStateSI& s)
        {
            const Int nelem = s.cfg.nelem;
            const Int npoin = s.cfg.npoin;
            const Int nep = s.cfg.nep;

            // node -> incident elements (CSR)
            std::vector<Int> cnt(static_cast<std::size_t>(npoin) + 2, 0);
            for (Int ie = 1; ie <= nelem; ++ie)
                for (Int in = 1; in <= nep; ++in)
                    ++cnt[static_cast<std::size_t>(s.intma(in, ie))];
            std::vector<Int> n2e_ptr(static_cast<std::size_t>(npoin) + 2, 0);
            for (Int ip = 1; ip <= npoin; ++ip)
                n2e_ptr[static_cast<std::size_t>(ip) + 1] =
                    n2e_ptr[static_cast<std::size_t>(ip)] + cnt[static_cast<std::size_t>(ip)];
            std::vector<Int> n2e(static_cast<std::size_t>(n2e_ptr[static_cast<std::size_t>(npoin) + 1]));
            std::vector<Int> fill(static_cast<std::size_t>(npoin) + 2, 0);
            for (Int ie = 1; ie <= nelem; ++ie)
                for (Int in = 1; in <= nep; ++in)
                {
                    const Int ip = s.intma(in, ie);
                    const std::size_t pos =
                        static_cast<std::size_t>(n2e_ptr[static_cast<std::size_t>(ip)] + fill[static_cast<std::size_t>(ip)]);
                    n2e[pos] = ie;
                    ++fill[static_cast<std::size_t>(ip)];
                }

            // greedy coloring
            std::vector<Int> color(static_cast<std::size_t>(nelem) + 1, -1);
            std::vector<Int> used;            // color-id marks, stamped per element
            std::vector<Int> stamp;           // forbidden[c]==ie when neighbour uses c
            Int ncolor = 0;

            for (Int ie = 1; ie <= nelem; ++ie)
            {
                // mark colors used by node-neighbours of ie
                for (Int in = 1; in <= nep; ++in)
                {
                    const Int ip = s.intma(in, ie);
                    for (Int k = n2e_ptr[static_cast<std::size_t>(ip)];
                         k < n2e_ptr[static_cast<std::size_t>(ip) + 1]; ++k)
                    {
                        const Int je = n2e[static_cast<std::size_t>(k)];
                        const Int cj = color[static_cast<std::size_t>(je)];
                        if (cj >= 0)
                        {
                            if (cj >= static_cast<Int>(stamp.size())) stamp.resize(cj + 1, -1);
                            stamp[static_cast<std::size_t>(cj)] = ie;
                        }
                    }
                }
                Int c = 0;
                while (c < static_cast<Int>(stamp.size()) && stamp[static_cast<std::size_t>(c)] == ie) ++c;
                color[static_cast<std::size_t>(ie)] = c;
                if (c + 1 > ncolor) ncolor = c + 1;
            }

            // group elements by color (counting sort)
            s.ncolor = ncolor;
            s.color_ptr.assign(static_cast<std::size_t>(ncolor) + 1, 0);
            for (Int ie = 1; ie <= nelem; ++ie)
                ++s.color_ptr[static_cast<std::size_t>(color[static_cast<std::size_t>(ie)]) + 1];
            for (Int c = 0; c < ncolor; ++c)
                s.color_ptr[static_cast<std::size_t>(c) + 1] += s.color_ptr[static_cast<std::size_t>(c)];
            s.color_elem.assign(static_cast<std::size_t>(nelem), 0);
            std::vector<Int> off(s.color_ptr.begin(), s.color_ptr.end());
            for (Int ie = 1; ie <= nelem; ++ie)
            {
                const Int c = color[static_cast<std::size_t>(ie)];
                s.color_elem[static_cast<std::size_t>(off[static_cast<std::size_t>(c)]++)] = ie;
            }

            // Build the pressure-only coloured list used in the pressure CG
            // matrix-vector product.  CHT pressure/momentum are fluid-only, so
            // visiting solid tetrahedra in every pressure matvec is pure cost.
            // We keep the same colours as the full list; therefore the race-free
            // nodal scatter property is preserved without recolouring.
            s.pressure_ncolor = ncolor;
            s.pressure_color_ptr.assign(static_cast<std::size_t>(ncolor) + 1, 0);

            for (Int c = 0; c < ncolor; ++c)
            {
                const Int beg = s.color_ptr[static_cast<std::size_t>(c)];
                const Int end = s.color_ptr[static_cast<std::size_t>(c) + 1];

                for (Int k = beg; k < end; ++k)
                {
                    const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                    if (s.mat_elem(ie) == 0)
                    {
                        ++s.pressure_color_ptr[static_cast<std::size_t>(c) + 1];
                    }
                }
            }

            for (Int c = 0; c < ncolor; ++c)
            {
                s.pressure_color_ptr[static_cast<std::size_t>(c) + 1] +=
                    s.pressure_color_ptr[static_cast<std::size_t>(c)];
            }

            s.pressure_nelem = s.pressure_color_ptr[static_cast<std::size_t>(ncolor)];
            s.pressure_color_elem.assign(static_cast<std::size_t>(s.pressure_nelem), 0);

            std::vector<Int> pressure_off(s.pressure_color_ptr.begin(), s.pressure_color_ptr.end());

            for (Int c = 0; c < ncolor; ++c)
            {
                const Int beg = s.color_ptr[static_cast<std::size_t>(c)];
                const Int end = s.color_ptr[static_cast<std::size_t>(c) + 1];

                for (Int k = beg; k < end; ++k)
                {
                    const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                    if (s.mat_elem(ie) == 0)
                    {
                        const std::size_t pos =
                            static_cast<std::size_t>(pressure_off[static_cast<std::size_t>(c)]++);
                        s.pressure_color_elem[pos] = ie;
                    }
                }
            }

            (void)used;
        }
    };
}
