/**
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "PythonBindings/HHOInterface.h"

#include "aster_pybind.h"

void exportHHOToPython( py::module_ &mod ) {

    py::class_< HHO, HHO::HHOPtr >( mod, "HHO" )
        .def( py::init( &initFactoryPtr< HHO, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "projectOnLagrangeSpace", &HHO::projectOnLagrangeSpace,
              R"(
      Project field from HHO-space to Lagrange-space

      Arguments:
            hho_field (FieldOnNodesReal): hho field like displacement or thermic

      Returns:
            FieldOnNodesReal: HHO field project on Lagrange space
        )",
              py::arg( "hho_field" ) )
        .def( "projectOnHHOSpace", &HHO::projectOnHHOSpace,
              R"(
      Project real value to HHO-space

      Arguments:
            value (float): value to project

      Returns:
            FieldOnNodesReal: HHO field
        )",
              py::arg( "value" ) );
};
