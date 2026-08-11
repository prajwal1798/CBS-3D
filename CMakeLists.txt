cmake_minimum_required(VERSION 3.16)

project(CBS3Dpp_SI LANGUAGES CXX)

#==============================================================================
# Optional solver capabilities
#==============================================================================
option(CBS3D_ENABLE_OPENMP "Enable OpenMP shared-memory parallelism" ON)
option(CBS3D_ENABLE_MPI "Build the MPI distributed application" OFF)
option(CBS3D_ENABLE_PETSC "Enable the PETSc pressure-solver path" OFF)
option(
    CBS3D_KEEP_LEGACY_EXECUTABLE_NAME
    "Write the selected application as cbs3dpp_si for existing run scripts"
    ON
)
option(
    CBS3D_BUILD_SOURCE_VALIDATION
    "Build the Step-4 solid volumetric-source validation"
    ON
)

if (CBS3D_ENABLE_MPI AND NOT CBS3D_ENABLE_PETSC)
    message(FATAL_ERROR
        "The production MPI application requires PETSc. Configure with "
        "CBS3D_ENABLE_MPI=ON and CBS3D_ENABLE_PETSC=ON.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# Shared numerical source tree
#==============================================================================
set(CBS3D_CORE_SOURCES
    src/io/MeshIO.cpp
    src/io/Post.cpp
    src/io/DistributedPost.cpp
    src/io/RestartIO.cpp

    src/preprocess/Preprocess.cpp
    src/boundary/Boundary.cpp
    src/timestep/TimeStep.cpp

    src/solver/Convergence.cpp
    src/solver/Steps.cpp
    src/solver/Solver.cpp
    src/solver/SolverDistributedStep2.cpp
    src/solver/SolverDistributedStep3.cpp
    src/solver/SolverDistributedLoop.cpp
    src/solver/SolverDistributedProductionLoopRestartAfc.cpp
    src/solver/SolverDistributedWallAudit.cpp
    src/solver/SolverDistributedPressureAudit.cpp
    src/solver/SolverDistributedVelocityBoundaryInventory.cpp

    src/linalg/BandedMatrix.cpp
    src/linalg/BandedGaussianSolver.cpp
    src/linalg/CSRMatrix.cpp
    src/linalg/MatrixVectorCalc.cpp
    src/linalg/ConjugateGradient.cpp

    src/assembly/PressureAssembly.cpp
    src/assembly/VelocityCorrectionAssembly.cpp
    src/assembly/MomentumAssembly.cpp
    src/assembly/EnergyAssembly.cpp
    src/assembly/TemperatureAFC.cpp
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
    list(APPEND CBS3D_CORE_SOURCES
        src/linalg/PetscPressureSolver.cpp
        src/linalg/PetscDistributedPressureSolver.cpp
        src/linalg/PetscPersistentDistributedPressureSystem.cpp
    )
endif()

#==============================================================================
# Application selection
#==============================================================================
if (CBS3D_ENABLE_MPI)
    set(CBS3D_APPLICATION_TARGET cbs3d_parallel)
    set(CBS3D_APPLICATION_SOURCE apps/parallel/main.cpp)
else()
    set(CBS3D_APPLICATION_TARGET cbs3d_serial)
    set(CBS3D_APPLICATION_SOURCE apps/serial/main.cpp)
endif()

add_executable(
    ${CBS3D_APPLICATION_TARGET}
    ${CBS3D_APPLICATION_SOURCE}
    ${CBS3D_CORE_SOURCES}
)

if (CBS3D_KEEP_LEGACY_EXECUTABLE_NAME)
    set_target_properties(
        ${CBS3D_APPLICATION_TARGET}
        PROPERTIES OUTPUT_NAME cbs3dpp_si
    )
endif()

target_include_directories(${CBS3D_APPLICATION_TARGET} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_compile_features(${CBS3D_APPLICATION_TARGET} PRIVATE cxx_std_20)

#==============================================================================
# Compiler warnings and optimisation
#==============================================================================
if (MSVC)
    target_compile_options(${CBS3D_APPLICATION_TARGET} PRIVATE
        /W4
        /permissive-
    )
else()
    target_compile_options(${CBS3D_APPLICATION_TARGET} PRIVATE
        -O2
        -Wall
        -Wextra
        -pedantic
    )
endif()

#==============================================================================
# Step-4 solid volumetric-source validation
#==============================================================================
if (CBS3D_BUILD_SOURCE_VALIDATION)
    enable_testing()

    add_executable(cbs3d_source_validation
        tests/Step4SourceSolidValidation.cpp
        src/assembly/EnergyAssembly.cpp
    )

    target_include_directories(cbs3d_source_validation PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )

    target_compile_features(cbs3d_source_validation PRIVATE cxx_std_20)

    if (MSVC)
        target_compile_options(cbs3d_source_validation PRIVATE
            /W4
            /permissive-
        )
    else()
        target_compile_options(cbs3d_source_validation PRIVATE
            -O2
            -Wall
            -Wextra
            -pedantic
        )
    endif()

    add_test(
        NAME step4_solid_volumetric_source
        COMMAND cbs3d_source_validation
    )
endif()

#==============================================================================
# OpenMP
#==============================================================================
if (CBS3D_ENABLE_OPENMP)
    find_package(OpenMP)

    if (OpenMP_CXX_FOUND)
        target_link_libraries(
            ${CBS3D_APPLICATION_TARGET}
            PRIVATE OpenMP::OpenMP_CXX
        )

        target_compile_definitions(
            ${CBS3D_APPLICATION_TARGET}
            PRIVATE CBS3D_USE_OPENMP=1
        )
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

    target_link_libraries(
        ${CBS3D_APPLICATION_TARGET}
        PRIVATE MPI::MPI_CXX
    )

    target_compile_definitions(
        ${CBS3D_APPLICATION_TARGET}
        PRIVATE CBS3D_USE_MPI=1
    )
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

    target_compile_definitions(
        ${CBS3D_APPLICATION_TARGET}
        PRIVATE CBS3D_USE_PETSC=1
    )

    target_include_directories(${CBS3D_APPLICATION_TARGET} PRIVATE
        "${PETSC_DIR}/include"
        "${PETSC_PREFIX}/include"
    )

    target_link_libraries(${CBS3D_APPLICATION_TARGET} PRIVATE
        "${PETSC_LIBRARY}"
    )

    if (UNIX AND NOT APPLE)
        set_target_properties(
            ${CBS3D_APPLICATION_TARGET}
            PROPERTIES BUILD_RPATH "${PETSC_PREFIX}/lib"
        )
    endif()
endif()

install(
    TARGETS ${CBS3D_APPLICATION_TARGET}
    RUNTIME DESTINATION bin
)

#==============================================================================
# Configuration summary
#==============================================================================
get_target_property(
    CBS3D_CONFIGURED_OUTPUT_NAME
    ${CBS3D_APPLICATION_TARGET}
    OUTPUT_NAME
)

if (NOT CBS3D_CONFIGURED_OUTPUT_NAME)
    set(CBS3D_CONFIGURED_OUTPUT_NAME ${CBS3D_APPLICATION_TARGET})
endif()

message(STATUS "CBS3D++_SI configuration")
message(STATUS "  Application target : ${CBS3D_APPLICATION_TARGET}")
message(STATUS "  Executable name     : ${CBS3D_CONFIGURED_OUTPUT_NAME}")
message(STATUS "  OpenMP              : ${CBS3D_ENABLE_OPENMP}")
message(STATUS "  MPI                 : ${CBS3D_ENABLE_MPI}")
message(STATUS "  PETSc               : ${CBS3D_ENABLE_PETSC}")
message(STATUS "  Source validation   : ${CBS3D_BUILD_SOURCE_VALIDATION}")
