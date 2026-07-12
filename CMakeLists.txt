cmake_minimum_required(VERSION 3.16)

project(CBS3Dpp_SI LANGUAGES CXX)

#==============================================================================
# Optional solver capabilities
#==============================================================================
option(CBS3D_ENABLE_OPENMP "Enable OpenMP shared-memory parallelism" ON)
option(CBS3D_ENABLE_MPI "Enable MPI partitioned-solver development" OFF)
option(CBS3D_ENABLE_PETSC "Enable the PETSc pressure-solver path" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# Source files
#==============================================================================
set(CBS3D_SOURCES
    src/main.cpp

    src/io/MeshIO.cpp
    src/io/Post.cpp

    src/preprocess/Preprocess.cpp
    src/boundary/Boundary.cpp
    src/timestep/TimeStep.cpp

    src/solver/Convergence.cpp
    src/solver/Steps.cpp
    src/solver/Solver.cpp

    src/linalg/BandedMatrix.cpp
    src/linalg/BandedGaussianSolver.cpp
    src/linalg/CSRMatrix.cpp
    src/linalg/MatrixVectorCalc.cpp
    src/linalg/ConjugateGradient.cpp

    src/assembly/PressureAssembly.cpp
    src/assembly/MomentumAssembly.cpp
    src/assembly/EnergyAssembly.cpp

    src/assembly/SpalartAllmarasAssembly.cpp

    src/turbulence/SpalartAllmaras.cpp
    src/turbulence/WallDistance.cpp
    src/turbulence/TurbulencePreprocess.cpp

    src/boundary/TurbulenceBoundary.cpp

    src/utils/SolverProfiler.cpp

    src/parallel/PartitionMetadata.cpp
    src/parallel/HaloExchange.cpp
)

if (CBS3D_ENABLE_PETSC)
    list(APPEND CBS3D_SOURCES
        src/linalg/PetscPressureSolver.cpp
    )
endif()

add_executable(cbs3dpp_si ${CBS3D_SOURCES})

target_include_directories(cbs3dpp_si PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_compile_features(cbs3dpp_si PRIVATE cxx_std_20)

#==============================================================================
# Compiler warnings and optimisation
#==============================================================================
if (MSVC)
    target_compile_options(cbs3dpp_si PRIVATE
        /W4
        /permissive-
    )
else()
    target_compile_options(cbs3dpp_si PRIVATE
        -O2
        -Wall
        -Wextra
        -pedantic
    )
endif()

#==============================================================================
# OpenMP
#==============================================================================
if (CBS3D_ENABLE_OPENMP)
    find_package(OpenMP)

    if (OpenMP_CXX_FOUND)
        target_link_libraries(cbs3dpp_si PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(cbs3dpp_si PRIVATE CBS3D_USE_OPENMP=1)
    else()
        message(WARNING
            "OpenMP was requested but no C++ OpenMP implementation was found. "
            "Building without OpenMP.")
    endif()
endif()

#==============================================================================
# MPI
#==============================================================================
if (CBS3D_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)

    target_link_libraries(cbs3dpp_si PRIVATE MPI::MPI_CXX)
    target_compile_definitions(cbs3dpp_si PRIVATE CBS3D_USE_MPI=1)
endif()

#==============================================================================
# PETSc
#==============================================================================
if (CBS3D_ENABLE_PETSC)
    if (NOT DEFINED ENV{PETSC_DIR})
        message(FATAL_ERROR
            "CBS3D_ENABLE_PETSC=ON requires PETSC_DIR to point to PETSc.")
    endif()

    set(PETSC_DIR "$ENV{PETSC_DIR}")
    set(PETSC_PREFIX "${PETSC_DIR}")

    if (DEFINED ENV{PETSC_ARCH}
        AND EXISTS "${PETSC_DIR}/$ENV{PETSC_ARCH}/include")
        set(PETSC_PREFIX "${PETSC_DIR}/$ENV{PETSC_ARCH}")
    endif()

    find_library(PETSC_LIBRARY
        NAMES petsc
        HINTS
            "${PETSC_PREFIX}/lib"
            "${PETSC_DIR}/lib"
        REQUIRED
    )

    target_compile_definitions(cbs3dpp_si PRIVATE CBS3D_USE_PETSC=1)

    target_include_directories(cbs3dpp_si PRIVATE
        "${PETSC_DIR}/include"
        "${PETSC_PREFIX}/include"
    )

    target_link_libraries(cbs3dpp_si PRIVATE
        "${PETSC_LIBRARY}"
    )

    if (UNIX AND NOT APPLE)
        set_target_properties(cbs3dpp_si PROPERTIES
            BUILD_RPATH "${PETSC_PREFIX}/lib"
        )
    endif()
endif()

#==============================================================================
# Configuration summary
#==============================================================================
message(STATUS "CBS3D++_SI configuration")
message(STATUS "  OpenMP : ${CBS3D_ENABLE_OPENMP}")
message(STATUS "  MPI    : ${CBS3D_ENABLE_MPI}")
message(STATUS "  PETSc  : ${CBS3D_ENABLE_PETSC}")
