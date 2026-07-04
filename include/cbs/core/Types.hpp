#pragma once

//=============================================================================
// CBS3D++_SI
//
// Fundamental scalar and index types used throughout the solver.
//
// Int
//     Signed integer used for mesh indices, element numbers, node numbers,
//     boundary identifiers and iteration counters.
//
// Real
//     Floating-point type used for coordinates, physical variables,
//     finite-element coefficients and residuals.
//
// Size
//     Unsigned container-size type used when indexing std::vector storage.
//
// Mesh and finite-element arrays use one-based scientific indexing, while
// std::vector and other standard-library containers remain zero-based.
//=============================================================================

#include <cstddef>
#include <cstdint>

namespace cbs
{
    using Int = std::int32_t;
    using Real = double;
    using Size = std::size_t;
}
