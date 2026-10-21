#pragma once

//=============================================================================
// CBS3D++_SI
//
// Single entry point for PETSc headers.
//
// Some packaged PETSc builds are compiled without C++ support and therefore
// omit PETSC_FUNCTION_NAME_CXX from petscconf.h.  The public PETSc headers
// still select that macro whenever they are included from C++, which makes a
// plain #include <petscksp.h> fail with
//
//     'PETSC_FUNCTION_NAME_CXX' was not declared in this scope
//
// Define the standard value before the first PETSc header.  A PETSc build that
// already provides the macro keeps its own definition.
//=============================================================================

#if !defined(PETSC_FUNCTION_NAME_CXX)
#define PETSC_FUNCTION_NAME_CXX __func__
#endif

#if !defined(PETSC_CXX_RESTRICT)
#if defined(_MSC_VER)
#define PETSC_CXX_RESTRICT __restrict
#else
#define PETSC_CXX_RESTRICT __restrict__
#endif
#endif

#include <petscksp.h>
