/**
 * @file MedCouplingConversionInterface.cxx
 * @brief Python bindings for MedCouplingConversion.
 * @author Francesco Bettonte
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: francesco.bettonte@edf.fr */

#include "aster_pybind.h"

#include <MedUtils/MedCouplingConversion.h>

void exportMedCouplingConversionToPython( py::module_ &mod ) {

    mod.def( "getMedCouplingConversionData", &getMedCouplingConversionData, R"(
Return three dictionnaries containing data to create an equivalent MedCoupling unstructured mesh.

MedCoupling needs a mesh splitted by dimension for what concerns cells and groups of cells.
The group of nodes all belongs to an unique level so there is no need to split them.

 - The first dictionnary (cells) contains for each dimension (the keys) :
   1. The connectivity
   2. The connectivity index
 - The second dictionnary (groups_c) contains for each dimension (the keys) a dictionnary
   which keys are the groups names at the items their cells.
 - The third dictionnary (groups_n) contains for each group of nodes (the keys)
   the nodes composing the group.

Arguments:
    mesh (BaseMeshPtr): The aster mesh to be processed.

Returns:
    tuple (cells, groups_c, groups_n) : The data to create the equivalent MedCoupling mesh.
        )",
             py::arg( "mesh" ) );
}
