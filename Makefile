#==============================================================================
# CBS3D++_SI Makefile
#
# Supported build modes:
#
#   Serial/OpenMP:
#       make
#
#   MPI/OpenMP bootstrap:
#       make USE_MPI=1
#
#   Serial/OpenMP with PETSc:
#       make USE_PETSC=1 PETSC_DIR=/path/to/petsc
#
#   MPI/OpenMP with PETSc:
#       make USE_MPI=1 USE_PETSC=1 PETSC_DIR=/path/to/petsc
#
# Run examples:
#
#       make run CASE=simple
#
#       make run-mpi USE_MPI=1 NP=4 CASE=simple
#
# On Sunbird, load the required compiler and MPI modules before building:
#
#       module load <gcc-module>
#       module load openmpi/4.1.6
#
# The MPI compiler wrapper may be changed with:
#
#       make USE_MPI=1 MPI_CXX=mpicxx
#
# Disable OpenMP when required:
#
#       make USE_OPENMP=0
#
# Build variants use separate object directories so that serial, MPI, PETSc and
# OpenMP flags cannot accidentally reuse incompatible object files.
#==============================================================================

SERIAL_CXX ?= g++
MPI_CXX    ?= mpicxx
MPIEXEC    ?= mpirun

CASE           ?= simple
NP             ?= 2
PARTITION_ROOT ?= cbs_partitions_$(NP)

USE_OPENMP ?= 1
USE_MPI    ?= 0
USE_PETSC  ?= 0

#------------------------------------------------------------------------------
# Compiler selection
#------------------------------------------------------------------------------
ifeq ($(USE_MPI),1)
    CXX := $(MPI_CXX)
else
    CXX := $(SERIAL_CXX)
endif

#------------------------------------------------------------------------------
# Build-variant names
#------------------------------------------------------------------------------
BUILD_TAG := serial
TARGET    := cbs3dpp_si

ifeq ($(USE_MPI),1)
    BUILD_TAG := mpi
    TARGET    := $(TARGET)_mpi
endif

ifeq ($(USE_PETSC),1)
    BUILD_TAG := $(BUILD_TAG)-petsc
    TARGET    := $(TARGET)_petsc
endif

ifeq ($(USE_OPENMP),1)
    BUILD_TAG := $(BUILD_TAG)-openmp
else
    BUILD_TAG := $(BUILD_TAG)-noopenmp
    TARGET    := $(TARGET)_noopenmp
endif

BUILD_DIR := build/$(BUILD_TAG)
OBJ_DIR   := $(BUILD_DIR)/obj

#------------------------------------------------------------------------------
# Compilation and linking flags
#------------------------------------------------------------------------------
CPPFLAGS := -Iinclude
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -pedantic
LDFLAGS  ?=
LDLIBS   ?=

ifeq ($(USE_OPENMP),1)
    CXXFLAGS += -fopenmp
    LDFLAGS  += -fopenmp
    CPPFLAGS += -DCBS3D_USE_OPENMP=1
endif

ifeq ($(USE_MPI),1)
    CPPFLAGS += -DCBS3D_USE_MPI=1
endif

ifeq ($(USE_PETSC),1)
    ifndef PETSC_DIR
        $(error USE_PETSC=1 requires PETSC_DIR to point to the PETSc installation)
    endif

    PETSC_PREFIX := $(PETSC_DIR)

    ifdef PETSC_ARCH
        ifneq ($(wildcard $(PETSC_DIR)/$(PETSC_ARCH)/include),)
            PETSC_PREFIX := $(PETSC_DIR)/$(PETSC_ARCH)
        endif
    endif

    CPPFLAGS += -DCBS3D_USE_PETSC=1
    CPPFLAGS += -I$(PETSC_DIR)/include
    CPPFLAGS += -I$(PETSC_PREFIX)/include

    LDFLAGS += -L$(PETSC_PREFIX)/lib
    LDFLAGS += -Wl,-rpath,$(PETSC_PREFIX)/lib

    LDLIBS += -lpetsc
endif

#------------------------------------------------------------------------------
# Source files
#------------------------------------------------------------------------------
SRC := \
    src/main.cpp \
    src/io/MeshIO.cpp \
    src/io/Post.cpp \
    src/preprocess/Preprocess.cpp \
    src/boundary/Boundary.cpp \
    src/timestep/TimeStep.cpp \
    src/solver/Convergence.cpp \
    src/solver/Steps.cpp \
    src/solver/Solver.cpp \
    src/linalg/BandedMatrix.cpp \
    src/linalg/BandedGaussianSolver.cpp \
    src/linalg/CSRMatrix.cpp \
    src/linalg/MatrixVectorCalc.cpp \
    src/linalg/ConjugateGradient.cpp \
    src/assembly/PressureAssembly.cpp \
    src/assembly/MomentumAssembly.cpp \
    src/assembly/EnergyAssembly.cpp \
    src/utils/SolverProfiler.cpp \
    src/parallel/PartitionMetadata.cpp \
    src/parallel/HaloExchange.cpp

ifeq ($(USE_PETSC),1)
    SRC += src/linalg/PetscPressureSolver.cpp
endif

OBJ := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC))
DEP := $(OBJ:.o=.d)

#------------------------------------------------------------------------------
# Build rules
#------------------------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $@
	@echo
	@echo "Built $(TARGET)"
	@echo "  compiler : $(CXX)"
	@echo "  variant  : $(BUILD_TAG)"
	@echo

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

#------------------------------------------------------------------------------
# Run helpers
#------------------------------------------------------------------------------
run: $(TARGET)
	./$(TARGET) $(CASE)

ifeq ($(USE_MPI),1)
run-mpi: $(TARGET)
	$(MPIEXEC) -n $(NP) ./$(TARGET) $(CASE) $(PARTITION_ROOT)
else
run-mpi:
	@echo "run-mpi requires USE_MPI=1"
	@echo "Example: make run-mpi USE_MPI=1 NP=4 CASE=$(CASE)"
	@exit 1
endif

print-config:
	@echo "TARGET       = $(TARGET)"
	@echo "BUILD_TAG    = $(BUILD_TAG)"
	@echo "CXX          = $(CXX)"
	@echo "USE_OPENMP   = $(USE_OPENMP)"
	@echo "USE_MPI      = $(USE_MPI)"
	@echo "USE_PETSC    = $(USE_PETSC)"
	@echo "PETSC_DIR    = $(PETSC_DIR)"
	@echo "PETSC_ARCH   = $(PETSC_ARCH)"
	@echo "PARTITION_ROOT = $(PARTITION_ROOT)"

clean:
	rm -rf build
	rm -f cbs3dpp_si \
	      cbs3dpp_si_mpi \
	      cbs3dpp_si_petsc \
	      cbs3dpp_si_mpi_petsc \
	      cbs3dpp_si_noopenmp \
	      cbs3dpp_si_mpi_noopenmp \
	      cbs3dpp_si_petsc_noopenmp \
	      cbs3dpp_si_mpi_petsc_noopenmp

-include $(DEP)

.PHONY: all run run-mpi print-config clean
