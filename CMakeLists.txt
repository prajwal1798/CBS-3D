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
option(CBS3D_BUILD_TESTS "Build the verification and regression tests" OFF)
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
    # StepsWallModel.cpp compiles the validated Steps.cpp implementation through
    # narrow wrappers.  With CBS3D_SA_WALL_TREATMENT unset, those wrappers are
    # exact no-ops and retain the established wall-resolved behaviour.
    src/solver/StepsWallModel.cpp
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

    # RANS/CBS Step 1 and SA transport are taken directly from the
    # Liu/Nithiarasu P1-TET4 matrices in Chun-Bin Liu (2005), Appendix B.
    # Keep the historical experimental assemblies in the tree for A/B audit
    # work, but do not compile them into the production solver.
    src/assembly/MomentumAssemblyLiu.cpp
    src/assembly/SpalartAllmarasAssemblyLiu.cpp

    src/assembly/EnergyAssembly.cpp
    src/assembly/TemperatureAFC.cpp

    src/turbulence/SpalartAllmaras.cpp
    src/turbulence/WallDistance.cpp
    src/turbulence/TurbulencePreprocess.cpp
    src/turbulence/WallTreatment.cpp
    src/turbulence/WallModelCoupling.cpp
    src/turbulence/CHTWallModelCoupling.cpp

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
#
# MSVC's legacy /openmp mode implements only OpenMP 2.0 and rejects the
# atomic-write pragmas used by the coloured CSM assemblies.  `/openmp:llvm`
# provides the OpenMP 3.1+ runtime (libomp) that those kernels require.
#==============================================================================
function(cbs3d_link_openmp target)
    if (NOT CBS3D_ENABLE_OPENMP)
        return()
    endif()

    if (MSVC)
        target_compile_options(${target} PRIVATE /openmp:llvm)
        target_compile_definitions(${target} PRIVATE CBS3D_USE_OPENMP=1)
        target_link_libraries(${target} PRIVATE libomp)
        set(CBS3D_OPENMP_CXX_READY TRUE PARENT_SCOPE)
    else()
        find_package(OpenMP QUIET)
        if (OpenMP_CXX_FOUND)
            target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
            target_compile_definitions(${target} PRIVATE CBS3D_USE_OPENMP=1)
            set(CBS3D_OPENMP_CXX_READY TRUE PARENT_SCOPE)
        else()
            message(WARNING
                "OpenMP was requested but no C++ OpenMP implementation was found. "
                "Building without OpenMP.")
            set(CBS3D_OPENMP_CXX_READY FALSE PARENT_SCOPE)
        endif()
    endif()
endfunction()

cbs3d_link_openmp(${CBS3D_APPLICATION_TARGET})

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
    set(PETSC_LIB_HINTS "${PETSC_DIR}/lib")
    set(PETSC_INCLUDE_DIRS
        "${PETSC_DIR}/include"
        "${PETSC_DIR}/include/petsc"
    )

    if (DEFINED ENV{PETSC_ARCH})
        set(PETSC_ARCH "$ENV{PETSC_ARCH}")

        # Source-tree or legacy install: ${PETSC_DIR}/${PETSC_ARCH}/lib.
        if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include")
            set(PETSC_PREFIX "${PETSC_DIR}/${PETSC_ARCH}")
            list(PREPEND PETSC_LIB_HINTS "${PETSC_PREFIX}/lib")
            list(PREPEND PETSC_INCLUDE_DIRS "${PETSC_PREFIX}/include")
        endif()

        # Packaged install (for example MSYS2 UCRT64): architecture headers
        # live in include/petsc/${PETSC_ARCH} and the import library in
        # lib/petsc/${PETSC_ARCH}.
        if (EXISTS "${PETSC_DIR}/include/petsc/${PETSC_ARCH}")
            list(PREPEND PETSC_INCLUDE_DIRS
                "${PETSC_DIR}/include/petsc/${PETSC_ARCH}"
            )
        endif()
        if (EXISTS "${PETSC_DIR}/lib/petsc/${PETSC_ARCH}")
            list(PREPEND PETSC_LIB_HINTS
                "${PETSC_DIR}/lib/petsc/${PETSC_ARCH}"
            )
            set(PETSC_PREFIX "${PETSC_DIR}/lib/petsc/${PETSC_ARCH}")
        endif()
    endif()

    find_library(PETSC_LIBRARY
        NAMES petsc
        HINTS ${PETSC_LIB_HINTS}
        REQUIRED
    )

    message(STATUS "  PETSc library       : ${PETSC_LIBRARY}")
    message(STATUS "  PETSc include dirs  : ${PETSC_INCLUDE_DIRS}")

    target_compile_definitions(
        ${CBS3D_APPLICATION_TARGET}
        PRIVATE CBS3D_USE_PETSC=1
    )

    target_include_directories(${CBS3D_APPLICATION_TARGET} PRIVATE
        ${PETSC_INCLUDE_DIRS}
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

#==============================================================================
# Verification and regression tests
#
# One switch builds every test that has a well-defined serial target:
#
#     cmake -S . -B build/tests -DCBS3D_BUILD_TESTS=ON
#     cmake --build build/tests
#     ctest --test-dir build/tests --output-on-failure
#
# The MPI-only regressions require a real MPI implementation (two ranks).
# They are opt-in because the packaged Windows PETSc used for the local tower
# is an MPIUNI build.
#==============================================================================
option(CBS3D_BUILD_MPI_TESTS "Build the two-rank MPI wall/topology regressions" OFF)

if (CBS3D_BUILD_TESTS)
    enable_testing()

    function(cbs3d_add_verification target)
        cmake_parse_arguments(V "MPI" "TEST_NAME" "SOURCES;ARGS" ${ARGN})

        if (V_MPI AND NOT (CBS3D_ENABLE_MPI AND TARGET MPI::MPI_CXX))
            return()
        endif()

        add_executable(${target} ${V_SOURCES})

        target_include_directories(${target} PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
        )

        target_compile_features(${target} PRIVATE cxx_std_20)

        if (MSVC)
            target_compile_options(${target} PRIVATE /W4 /permissive-)
        else()
            target_compile_options(${target} PRIVATE -O2 -Wall -Wextra -pedantic)
        endif()

        cbs3d_link_openmp(${target})

        if (V_MPI)
            target_compile_definitions(${target} PRIVATE CBS3D_USE_MPI=1)
            target_link_libraries(${target} PRIVATE MPI::MPI_CXX)
        endif()

        if (NOT V_TEST_NAME)
            set(V_TEST_NAME ${target})
        endif()

        if (V_MPI AND MPIEXEC_EXECUTABLE)
            add_test(
                NAME ${V_TEST_NAME}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                        $<TARGET_FILE:${target}> ${V_ARGS}
            )
        else()
            add_test(
                NAME ${V_TEST_NAME}
                COMMAND ${target} ${V_ARGS}
            )
        endif()
    endfunction()

    #--------------------------------------------------------------------------
    # Turbulence and wall-model kernels
    #--------------------------------------------------------------------------
    cbs3d_add_verification(wall_distance_verification
        TEST_NAME wall_distance_exactness
        ARGS 60 40 20
        SOURCES
            tests/WallDistanceVerification.cpp
            src/turbulence/WallDistance.cpp
    )

    cbs3d_add_verification(wall_treatment_verification
        TEST_NAME wall_treatment
        SOURCES
            tests/WallTreatmentVerification.cpp
            src/turbulence/WallTreatment.cpp
    )

    cbs3d_add_verification(wall_treatment_nasa_reference
        TEST_NAME wall_treatment_nasa_reference
        SOURCES
            tests/WallTreatmentNASAReference.cpp
            src/turbulence/WallTreatment.cpp
    )

    cbs3d_add_verification(wall_treatment_reynolds_sweep
        TEST_NAME wall_treatment_reynolds_sweep
        SOURCES
            tests/WallTreatmentReynoldsSweep.cpp
            src/turbulence/WallTreatment.cpp
    )

    cbs3d_add_verification(thermal_wall_treatment_verification
        TEST_NAME thermal_wall_treatment
        SOURCES tests/ThermalWallTreatmentVerification.cpp
    )

    cbs3d_add_verification(turbulence_wall_topology_verification
        TEST_NAME turbulence_wall_topology
        SOURCES
            tests/TurbulenceWallTopologyVerification.cpp
            src/turbulence/WallDistance.cpp
    )

    cbs3d_add_verification(wall_model_coupling_verification
        TEST_NAME wall_model_coupling
        SOURCES
            tests/WallModelCouplingVerification.cpp
            src/turbulence/WallModelCoupling.cpp
            src/turbulence/WallTreatment.cpp
    )

    cbs3d_add_verification(cht_wall_model_coupling_verification
        TEST_NAME cht_wall_model_coupling
        SOURCES
            tests/CHTWallModelCouplingVerification.cpp
            src/turbulence/CHTWallModelCoupling.cpp
            src/turbulence/WallModelCoupling.cpp
            src/turbulence/WallTreatment.cpp
    )

    #--------------------------------------------------------------------------
    # Spalart-Allmaras finite-element verification
    #--------------------------------------------------------------------------
    cbs3d_add_verification(sa_neg_element_verification
        TEST_NAME sa_neg_negative_branch
        SOURCES
            tests/SANegElementVerification.cpp
            src/assembly/SpalartAllmarasAssemblyLiu.cpp
            src/turbulence/SpalartAllmaras.cpp
    )

    cbs3d_add_verification(liu_nithiarasu_element_verification
        TEST_NAME liu_nithiarasu_element
        SOURCES
            tests/LiuNithiarasuElementVerification.cpp
            src/assembly/MomentumAssemblyLiu.cpp
            src/assembly/SpalartAllmarasAssemblyLiu.cpp
            src/turbulence/SpalartAllmaras.cpp
    )

    #--------------------------------------------------------------------------
    # Two-rank MPI regressions (opt-in; see CBS3D_BUILD_MPI_TESTS)
    #--------------------------------------------------------------------------
    if (CBS3D_BUILD_MPI_TESTS)
        cbs3d_add_verification(turbulence_wall_topology_mpi
            MPI
            TEST_NAME turbulence_wall_topology_mpi
            SOURCES tests/TurbulenceWallTopologyMPIVerification.cpp
        )

        cbs3d_add_verification(wall_model_coupling_mpi
            MPI
            TEST_NAME wall_model_coupling_mpi
            SOURCES
                tests/WallModelCouplingMPIRankLocalVerification.cpp
                src/turbulence/WallModelCoupling.cpp
                src/turbulence/WallTreatment.cpp
        )

        cbs3d_add_verification(cht_wall_model_coupling_mpi
            MPI
            TEST_NAME cht_wall_model_coupling_mpi
            SOURCES
                tests/CHTWallModelCouplingMPIVerification.cpp
                src/turbulence/CHTWallModelCoupling.cpp
                src/turbulence/WallModelCoupling.cpp
                src/turbulence/WallTreatment.cpp
        )
    endif()
endif()


