/**
 * @file SimpleFieldOnCellsInterface.cxx
 * @brief Interface python de SimpleFieldOnCells
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/SimpleFieldOnCellsInterface.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( getvalues_overloads, getValues, 0, 1 )

void exportSimpleFieldOnCellsToPython() {
    py::class_< SimpleFieldOnCellsReal, SimpleFieldOnCellsRealPtr,
                py::bases< DataStructure > >( "SimpleFieldOnCellsReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< SimpleFieldOnCellsReal >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< SimpleFieldOnCellsReal, std::string >))
        .def( "getValue", &SimpleFieldOnCellsReal::getValue,
              py::return_value_policy< py::return_by_value >(), R"(
Returns the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt` sub-point.

Args:
    ima  (int): Index of cells.
    icmp (int): Index of component.
    ipt  (int): Index of point.
    ispt (int): Index of sub-point.

Returns:
    (float): Value of field at *ima*, of *icmp*, at *ipt*, at *ispt*;
             NaN if the position is not allocated.
        )", (py::arg( "self" ), py::arg("ima"), py::arg("icmp"), py::arg("ipt"), py::arg("ispt")))

        .def( "getValues", &SimpleFieldOnCellsReal::getValues, getvalues_overloads( R"(
Returns two numpy arrays with shape ( number_of_cells_with_components, number_of_components )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Args:
        copy (bool): If True copy the data, default: *False*

Returns:
    ndarray (float): Field values.
    ndarray (bool): Mask for the field values.
        )", (py::arg( "self" ), py::arg("copy" ))))

        .def( "getCellsWithComponents", &SimpleFieldOnCellsReal::getCellsWithComponents, R"(
Returns the list of cells where the field is defined.

Returns:
    tuple (int): Indexes of cells where the field is defined.
        )", (py::arg( "self" )) )
        .def( "getNumberOfComponents", &SimpleFieldOnCellsReal::getNumberOfComponents )
        .def( "getNameOfComponent", &SimpleFieldOnCellsReal::getNameOfComponent )
        .def( "getNameOfComponents", &SimpleFieldOnCellsReal::getNameOfComponents )
        .def( "getNumberOfCells", &SimpleFieldOnCellsReal::getNumberOfCells )
        .def( "getMaxNumberOfPoints", &SimpleFieldOnCellsReal::getMaxNumberOfPoints )
        .def( "getMaxNumberOfSubPoints", &SimpleFieldOnCellsReal::getMaxNumberOfSubPoints )
        .def( "getNumberOfPointsOfCell", &SimpleFieldOnCellsReal::getNumberOfPointsOfCell )
        .def( "getNumberOfSubPointsOfCell", &SimpleFieldOnCellsReal::getNumberOfSubPointsOfCell )
        .def( "getNumberOfComponentsForSubpointsOfCell",
              &SimpleFieldOnCellsReal::getNumberOfComponentsForSubpointsOfCell )
        .def( "getPhysicalQuantity", &SimpleFieldOnCellsReal::getPhysicalQuantity )
        .def( "getFieldLocation", &SimpleFieldOnCellsReal::getFieldLocation )
        .def( "updateValuePointers", &SimpleFieldOnCellsReal::updateValuePointers );
};
