/**
 * @file AssemblyMatrixInterface.cxx
 * @brief Interface python de AssemblyMatrix
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
#include "PythonBindings/AssemblyMatrixInterface.h"
#include <PythonBindings/factory.h>

void exportAssemblyMatrixToPython() {

    py::class_< AssemblyMatrixDisplacementReal, AssemblyMatrixDisplacementRealPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixDisplacementReal", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< AssemblyMatrixDisplacementReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixDisplacementReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixDisplacementReal, PhysicalProblemPtr > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixDisplacementReal::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixDisplacementReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixDisplacementReal::build )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" );

    py::class_< AssemblyMatrixDisplacementComplex, AssemblyMatrixDisplacementComplexPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixDisplacementComplex",
                                                   py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< AssemblyMatrixDisplacementComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixDisplacementComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixDisplacementComplex::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixDisplacementComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixDisplacementComplex::build )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixDisplacementComplex::transposeConjugate )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementComplex::getNumberOfElementaryMatrix );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureReal, AssemblyMatrixTemperatureRealPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixTemperatureReal", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor( &initFactoryPtr< AssemblyMatrixTemperatureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixTemperatureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixTemperatureReal::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixTemperatureReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixTemperatureReal::build )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixTemperatureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureComplex, AssemblyMatrixTemperatureComplexPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixTemperatureComplex", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< AssemblyMatrixTemperatureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixTemperatureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixTemperatureComplex::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixTemperatureComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixTemperatureComplex::build )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureComplex::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixTemperatureComplex::transposeConjugate );

    py::class_< AssemblyMatrixPressureReal, AssemblyMatrixPressureRealPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixPressureReal", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor( &initFactoryPtr< AssemblyMatrixPressureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< AssemblyMatrixPressureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixPressureReal::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixPressureReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixPressureReal::build )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixPressureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixPressureComplex, AssemblyMatrixPressureComplexPtr,
                py::bases< BaseAssemblyMatrix > >( "AssemblyMatrixPressureComplex", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor( &initFactoryPtr< AssemblyMatrixPressureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixPressureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix", &AssemblyMatrixPressureComplex::appendElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixPressureComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixPressureComplex::build )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureComplex::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixPressureComplex::transposeConjugate );
};
