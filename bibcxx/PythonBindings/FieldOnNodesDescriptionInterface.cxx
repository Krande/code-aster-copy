/**
 * @file FieldOnNodesDescriptionInterface.cxx
 * @brief Interface python de FieldOnNodesDescription
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/FieldOnNodesDescriptionInterface.h"

#include "aster_pybind.h"

void exportFieldOnNodesDescriptionToPython( py::module_ &mod ) {

    py::class_< FieldOnNodesDescription, FieldOnNodesDescription::FieldOnNodesDescriptionPtr,
                DataStructure >( mod, "FieldOnNodesDescription" )
        .def( py::init( &initFactoryPtr< FieldOnNodesDescription > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesDescription, BaseMeshPtr > ) )
        .def( py::init(
            &initFactoryPtr< FieldOnNodesDescription, std::string, BaseMeshPtr, std::string > ) )
        .def( "getMesh", &FieldOnNodesDescription::getMesh )
        .def( "updateValuePointers", &FieldOnNodesDescription::updateValuePointers )
        .def( "getNumberOfDofs", &FieldOnNodesDescription::getNumberOfDofs, R"(
            Return the number of DOFs

            Returns:
                int: number of DOFs
      )" )
        .def( "getNodesAndComponentsFromDOF",
              &FieldOnNodesDescription::getNodesAndComponentsFromDOF,
              R"(
            Return the list of node id and name of component for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index

            Returns:
                list[tuple[int, str]] : node id and name of component for each dofs
            )",
              py::arg( "local" ) = true )
        .def( "getNodeAndComponentFromDOF", &FieldOnNodesDescription::getNodeAndComponentFromDOF,
              R"(
            Return the node id and name of component for the dof

            Arguments:
                dof (int): dof index
                local (bool) = True: if True use local node index else use global index

            Returns:
                list[tuple[int, str]] : node id and name of component for each dofs
            )",
              py::arg( "dof" ), py::arg( "local" ) = true )
        .def( "getNodesAndComponentsNumberFromDOF",
              &FieldOnNodesDescription::getNodesAndComponentsNumberFromDOF,
              R"(
            Return the list of node id and component id for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index

            Returns:
                list[tuple[int, int]] : node id and component if for each dofs
            )",
              py::arg( "local" ) = true )
        .def( "getNodeAndComponentNumberFromDOF",
              &FieldOnNodesDescription::getNodeAndComponentNumberFromDOF,
              R"(
            Return the node id and component id for dof

            Arguments:
                dof (int) : DOF indec
                local (bool) = True: if True use local node index else use global index

            Returns:
                list[tuple[int, int]] : node id and component if for each dofs
            )",
              py::arg( "dof" ), py::arg( "local" ) = true )
        .def( "getDOFsFromNodesAndComponentsNumber",
              &FieldOnNodesDescription::getDOFsFromNodesAndComponentsNumber,
              R"(
            Return the dict of dofs with the pair (node id, name id) as keys

            Arguments:
                local (bool) = True: if True use local node index else use global index

            Returns:
                dict[int, str] : dofs id for each node id and component id
            )",
              py::arg( "local" ) = true )
        .def( "getDOFsFromNodesAndComponents",
              &FieldOnNodesDescription::getDOFsFromNodesAndComponents,
              R"(
           Return the dict of dofs with the pair (node id, component's name) as keys

            Arguments:
                local (bool) = True: if True use local node index else use global index

            Returns:
                dict[int, str] : dofs id for each node id and component's name
            )",
              py::arg( "local" ) = true )
        .def( "getComponents", &FieldOnNodesDescription::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" );
};
