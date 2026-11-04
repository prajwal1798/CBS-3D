#pragma once

//=============================================================================
// CBS3D++_SI
//
// FILE      : PointSampler.hpp
// PURPOSE   : TIME HISTORIES AT FIXED POINTS DURING AN MPI RUN.
//
//             FIELDS (8): T, u, v, w, p, omega_x, omega_y, omega_z
//               T, u, p   P1 INTERPOLATION IN THE CONTAINING TETRAHEDRON
//                         (BARYCENTRIC WEIGHTS)
//               omega     curl u OF THAT ELEMENT (CONSTANT IN A P1 TET)
//             EACH POINT IS EVALUATED BY THE ONE RANK THAT OWNS ITS ELEMENT
//             (LOWEST RANK ON A SHARED FACE) AND COLLECTED ON RANK 0.
//
//   PROBES    CBS3D_PROBE_FILE=probes.txt     LINES "name x y z" (# COMMENTS)
//             CBS3D_PROBE_DT=0.01             SAMPLING INTERVAL [s] (TRANSIENT)
//             CBS3D_PROBE_EVERY=100           ITERATIONS (STEADY RUNS)
//             -> output/<case>/<case>_probes.csv
//                iteration,time,heat_scale,<name>_T,<name>_u,...,<name>_wz
//
//   GRID      CBS3D_GRID_DIMS="nx ny nz"      CELL-CENTRED POINTS OVER THE
//             CBS3D_GRID_BOUNDS="x0 x1 y0 y1 z0 z1"   (DEFAULT: MESH BOX)
//             CBS3D_GRID_DT=0.1               SAMPLING INTERVAL [s]
//             CBS3D_GRID_EVERY=1000           ITERATIONS (STEADY RUNS)
//             -> output/<case>/fno_grid/
//                  grid.json                  DIMENSIONS, BOUNDS, FIELD ORDER
//                  mask.i8                    int8 PER POINT: -1 OUTSIDE THE
//                                             MESH, ELSE THE MATERIAL ID
//                  index.csv                  file,iteration,time,heat_scale
//                  sample_<iteration>.f32     float32, LITTLE-ENDIAN,
//                                             [field][k][j][i], i FASTEST;
//                                             0 OUTSIDE THE MESH
//
//             A RUN THAT STARTS AT t = 0 ALSO SAMPLES ITS INITIAL STATE. A
//             RESUMED RUN APPENDS TO THE CSV FILES.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <string>

namespace cbs
{
    class PointSampler
    {
    public:
        static void initialise(CBSStateSI& s, const std::string& case_name);

        // SAMPLES WHATEVER IS DUE AT THE CURRENT TIME / ITERATION.
        // COLLECTIVE: EVERY RANK MUST CALL IT.
        static void sample(CBSStateSI& s, const std::string& case_name, Int iteration);

        static constexpr int nfields = 8;
    };
}
