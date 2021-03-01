/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/DOFNumberingInterface.h"
#include <PythonBindings/factory.h>

void exportDOFNumberingToPython() {


    py::class_< DOFNumberingClass, DOFNumberingClass::DOFNumberingPtr,
                py::bases< BaseDOFNumberingClass > >( "DOFNumbering", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< DOFNumberingClass >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DOFNumberingClass, std::string >))
        .def( "getDescription", &DOFNumberingClass::getDescription )
// -------------------------------------------------------------------------------------------------
        .def( "useLagrangeMultipliers", &DOFNumberingClass::useLagrangeMultipliers, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )",
              ( py::arg( "self" ) )  )
// -------------------------------------------------------------------------------------------------
        .def( "useSingleLagrangeMultipliers", &DOFNumberingClass::useSingleLagrangeMultipliers, R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )",
              ( py::arg( "self" ) )   )
// -------------------------------------------------------------------------------------------------
        .def( "getNodeAssociatedToRow", static_cast<ASTERINTEGER (DOFNumberingClass::*)
                (const ASTERINTEGER, bool) const> (&DOFNumberingClass::getNodeAssociatedToRow),
                R"(
Returns the node index associated to a dof index.

- If the row is associated to a physical DOF, the *positive* id of the node
  handling this DOF is returned.

- If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the *negative* id of the node which is constrained by
  the multiplier is returned.

- If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
  implying several DOF, a blank string ' ' is returned (since no component can be
  identified).

Arguments:
    row (int): Index of the dof.

    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering).

Returns:
    int: index of the dof.
        )",
              ( py::arg( "self"), py::args( "row", "local") ) )
// -------------------------------------------------------------------------------------------------
        .def( "getNodeAssociatedToRow", static_cast<ASTERINTEGER (DOFNumberingClass::*)
                   (const ASTERINTEGER) const> (&DOFNumberingClass::getNodeAssociatedToRow),
                   R"(
Returns the node index associated to a dof index.

- If the row is associated to a physical DOF, the *positive* id of the node
  handling this DOF is returned.

- If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the *negative* id of the node which is constrained by the
  the multiplier is returned.

- If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
  implying several DOF, a blank string ' ' is returned (since no component can be
  identified).

Arguments:
    row (int): Index of the dof.

Returns:
    int: index of the dof.
        )",
              ( py::arg( "self"), py::args( "row", "local") ) )
// -------------------------------------------------------------------------------------------------
        .def( "getRowsAssociatedToPhysicalDofs", static_cast<VectorLong (DOFNumberingClass::*)
                   () const> (&DOFNumberingClass::getRowsAssociatedToPhysicalDofs), R"(
Returns the indexes of the physical dof.

Returns:
    [int]: indexes of the physical dof.
        )",
              ( py::arg( "self" ) )  )
// -------------------------------------------------------------------------------------------------
        .def( "getRowsAssociatedToPhysicalDofs", static_cast<VectorLong (DOFNumberingClass::*)
                   (bool) const> (&DOFNumberingClass::getRowsAssociatedToPhysicalDofs), R"(
Returns the indexes of the physical dof.

Arguments:
    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    [int]: indexes of the physical dof.
        )",
              ( py::arg( "self" ) , py::arg( "local" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "getRowAssociatedToNodeComponent", static_cast<ASTERINTEGER (DOFNumberingClass::*)
                (const ASTERINTEGER, const std::string) const>
                    (&DOFNumberingClass::getRowAssociatedToNodeComponent),
                    R"(
Returns the index of the dof associated to a node.

Arguments:
    node (int): Index of the node.
    component (str): name of the component

Returns:
    int: index of the dof.
        )",
              ( py::arg( "self" ), py::args( "node", "component" )  )  )
// -------------------------------------------------------------------------------------------------
        .def( "getRowAssociatedToNodeComponent", static_cast<ASTERINTEGER (DOFNumberingClass::*)
                   (const ASTERINTEGER, const std::string, bool) const>
                    (&DOFNumberingClass::getRowAssociatedToNodeComponent),
                    R"(
Returns the index of the dof associated to a node.

Arguments:
    node (int): Index of the node.
    component (str): name of the component
    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    int: index of the dof.
        )",
              ( py::arg( "self" ), py::args( "node", "component", "local" )  )  )
// -------------------------------------------------------------------------------------------------
        .def( "getRowsAssociatedToLagrangeMultipliers",static_cast<VectorLong (DOFNumberingClass::*)
                   () const> (&DOFNumberingClass::getRowsAssociatedToLagrangeMultipliers), R"(
Returns the indexes of the Lagrange multipliers dof.

Returns:
    int: indexes of the Lagrange multipliers dof.
        )",
              ( py::arg( "self" ) )  )
        .def( "getComponents", &DOFNumberingClass::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )",
              ( py::arg( "self" ) )  )
// -------------------------------------------------------------------------------------------------
        .def( "getRowsAssociatedToLagrangeMultipliers",static_cast<VectorLong (DOFNumberingClass::*)
                   (bool) const> (&DOFNumberingClass::getRowsAssociatedToLagrangeMultipliers), R"(
Returns the indexes of the Lagrange multipliers dof.

Arguments:
    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    [int]: indexes of the Lagrange multipliers dof.
        )",
              ( py::arg( "self" ), py::arg( "local" ) )  )
// -------------------------------------------------------------------------------------------------
        .def( "getComponents", &DOFNumberingClass::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )",
              ( py::arg( "self" ) )  )
// -------------------------------------------------------------------------------------------------
        .def( "getComponentAssociatedToRow", static_cast<std::string (DOFNumberingClass::*)
                   (const ASTERINTEGER) const> (&DOFNumberingClass::getComponentAssociatedToRow),
                   R"(
Returns the component name associated to a dof index.

- If the row is associated to a physical DOF, the name of the component is returned.

- If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the component which is constrained by the multiplier is
  returned.

- If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
  implying several DOF, a blank string ' ' is returned (since no component can be
  identified).

Arguments:
    row (int): Index of the dof (1-based index).

Returns:
    str: component name.
        )",
              ( py::arg( "self" ), py::arg( "row" ))  )
// ------------------------------------------------------------------------------------------------
        .def( "getComponentAssociatedToRow", static_cast<std::string (DOFNumberingClass::*)
                (const ASTERINTEGER, bool) const> (&DOFNumberingClass::getComponentAssociatedToRow),
                R"(
Returns the component name associated to a dof index.

- If the row is associated to a physical DOF, the name of the component is returned.

- If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the component which is constrained by the multiplier is
  returned.

- If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
  implying several DOF, a blank string ' ' is returned (since no component can be
  identified).

Arguments:
    row (int): Index of the dof (1-based index).

    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    str: component name.
        )",
              ( py::arg( "self" ), py::args( "row", "local"))  )
// ------------------------------------------------------------------------------------------------
        .def( "getComponentsAssociatedToNode", static_cast<VectorString (DOFNumberingClass::*)
                   (const ASTERINTEGER) const> (&DOFNumberingClass::getComponentsAssociatedToNode),
                    R"(
Returns the components name associated to a node index.

Arguments:
    node (int): Index of the node (1-based index).

Returns:
    str: component names.
        )",
              ( py::arg( "self" ), py::arg( "node" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "getComponentsAssociatedToNode", static_cast<VectorString (DOFNumberingClass::*)
            (const ASTERINTEGER, bool) const> (&DOFNumberingClass::getComponentsAssociatedToNode),
            R"(
Returns the components name associated to a node index.

Arguments:
    node (int): Index of the node (1-based index).

    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    str: component names.
        )",
              ( py::arg( "self" ), py::args( "node" , "local")) )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfDofs", static_cast<ASTERINTEGER (DOFNumberingClass::*) () const>
                    (&DOFNumberingClass::getNumberOfDofs), R"(
Returns the number of DOFs.

Returns:
    int: number of DOFs.
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfDofs", static_cast<ASTERINTEGER (DOFNumberingClass::*) (const bool) const>
                   (&DOFNumberingClass::getNumberOfDofs), R"(
Returns the number of DOFs.

Arguments:
    local=false (bool): not used (kept for compatibility with ParallelDOFNumbering)

Returns:
    int: number of DOFs.
        )",
              ( py::arg( "self" ), py::arg( "local"))  );

};
