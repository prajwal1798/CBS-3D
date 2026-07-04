#pragma once

//=============================================================================
// CBS3D++_SI
//
// Input routines for the three-dimensional CBS solver.
//
// The following case files are read:
//
//     .plt       tetrahedral mesh and boundary faces
//     .bco       boundary-condition definitions
//     .par       solver and physical parameters
//     .material  material number assigned to each element
//     .matprop   material properties
//
// The material files are required only when material properties are enabled.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/core/CaseFiles.hpp"

#include <string>

namespace cbs
{
    class MeshIO
    {
    public:
        static void readAll(
            const std::string& case_name,
            CBSStateSI& state);

    private:
        static void readSizes(
            const CaseFiles& files,
            CBSStateSI& state);

        static void readMeshFile(
            const CaseFiles& files,
            CBSStateSI& state);

        static void readBoundaryFile(
            const CaseFiles& files,
            CBSStateSI& state);

        static void readParameterFile(
            const CaseFiles& files,
            CBSStateSI& state);

        static void readMaterialFile(
            const CaseFiles& files,
            CBSStateSI& state);

        static void readMaterialPropertyFile(
            const CaseFiles& files,
            CBSStateSI& state);

        static void initialiseFields(
            CBSStateSI& state);
    };
}