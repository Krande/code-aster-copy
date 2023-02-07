/**
 * @file SimpleFieldOnNodesInterface.cxx
 * @brief Interface python de SimpleFieldOnNodes
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

#include "PythonBindings/SimpleFieldOnNodesInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataStructureInterface.h"

void exportSimpleFieldOnNodesToPython( py::module_ &mod ) {
    py::class_< SimpleFieldOnNodesReal, SimpleFieldOnNodesRealPtr, DataStructure >(
        mod, "SimpleFieldOnNodesReal" )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesReal > ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesReal, std::string > ) )
        .def( "getValue", &SimpleFieldOnNodesReal::getValue, py::return_value_policy::copy,
              R"(
Returns the value of the `icmp` component of the field on the `ino` node.

Arguments:
        ino (int): Index of node.
        icmp (int): Index of component.

Returns:
    (float): The field value. NaN is returned if the position is not allocated.
        )",
              py::arg( "ino" ), py::arg( "icmp" ) )

        .def( "getValues", &SimpleFieldOnNodesReal::getValues, R"(
Returns two numpy arrays with shape ( number_of_components, space_dimension )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Args:
        copy (bool): If True copy the data, default: *False*

Returns:
    ndarray (float): Field values.
    ndarray (bool): Mask for the field values.
        )",
              ( py::arg( "self" ), py::arg( "copy" ) = false ) )

        .def( "getNumberOfComponents", &SimpleFieldOnNodesReal::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesReal::getNumberOfNodes )
        .def( "getNameOfComponents", &SimpleFieldOnNodesReal::getNameOfComponents )
        .def( "getNameOfComponent", &SimpleFieldOnNodesReal::getNameOfComponent )
        .def( "getPhysicalQuantity", &SimpleFieldOnNodesReal::getPhysicalQuantity )
        .def( "updateValuePointers", &SimpleFieldOnNodesReal::updateValuePointers );
    py::class_< SimpleFieldOnNodesComplex, SimpleFieldOnNodesComplexPtr, DataStructure >(
        mod, "SimpleFieldOnNodesComplex" )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesComplex > ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesComplex, std::string > ) )

        .def( "getValue", &SimpleFieldOnNodesComplex::getValue, py::return_value_policy::copy,
              R"(
Returns the value of the `icmp` component of the field on the `ino` node.

Arguments:
    ino (int): Index of node.
    icmp (int): Index of component.

Returns:
    complex: The field value. NaN is returned if the position is not allocated.
        )",
              py::arg( "ino" ), py::arg( "icmp" ) )

        .def( "getValues", &SimpleFieldOnNodesComplex::getValues,
              R"(
Returns two numpy arrays with shape ( number_of_components, space_dimension )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Args:
        copy (bool): If True copy the data, default: *False*

Returns:
    ndarray (complex): Field values.
    ndarray (bool): Mask for the field values.
        )",
              ( py::arg( "self" ), py::arg( "copy" ) = false ) )

        .def( "getNumberOfComponents", &SimpleFieldOnNodesComplex::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesComplex::getNumberOfNodes )
        .def( "getNameOfComponents", &SimpleFieldOnNodesComplex::getNameOfComponents )
        .def( "getNameOfComponent", &SimpleFieldOnNodesComplex::getNameOfComponent )
        .def( "getPhysicalQuantity", &SimpleFieldOnNodesComplex::getPhysicalQuantity )
        .def( "updateValuePointers", &SimpleFieldOnNodesComplex::updateValuePointers );
};
