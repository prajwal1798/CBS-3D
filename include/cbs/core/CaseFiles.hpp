#pragma once

//=============================================================================
// CBS3D++_SI
//
// File names associated with one solver case.
//
// A case is identified by one common base name. The complete input-file names
// are obtained by appending the required extension:
//
//     <case>.par        numerical and physical controls
//     <case>.plt        tetrahedral mesh and boundary faces
//     <case>.bco        boundary-flag to solver-BC mapping
//     <case>.var        restart solution fields
//     <case>.material   material number assigned to each element
//     <case>.matprop    density, heat capacity, conductivity, viscosity
//                      and volumetric heat-source properties
//
// The structure contains file names only. Opening, parsing and validation are
// performed by MeshIO.
//=============================================================================

#include <string>

namespace cbs
{
    struct CaseFiles
    {
        std::string case_name;

        std::string par;
        std::string plt;
        std::string bco;
        std::string var;
        std::string material;
        std::string matprop;

        explicit CaseFiles(const std::string& name)
            : case_name(name),
              par(name + ".par"),
              plt(name + ".plt"),
              bco(name + ".bco"),
              var(name + ".var"),
              material(name + ".material"),
              matprop(name + ".matprop")
        {
        }
    };
}
