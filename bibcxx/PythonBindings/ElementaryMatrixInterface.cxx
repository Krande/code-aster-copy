/**
 * @file ElementaryMatrixInterface.cxx
 * @brief Interface python de ElementaryMatrix
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

#include "PythonBindings/ElementaryMatrixInterface.h"

#include "aster_pybind.h"

void exportElementaryMatrixToPython( py::module_ &mod ) {

    py::class_< BaseElementaryMatrix, BaseElementaryMatrixPtr, DataStructure >(
        mod, "BaseElementaryMatrix" )
        // fake initFactoryPtr: not buildable
        // fake initFactoryPtr: not buildable
        .def( "getMaterialField", &BaseElementaryMatrix::getMaterialField )
        .def( "getElementaryCharacteristics", &BaseElementaryMatrix::getElementaryCharacteristics )
        .def( "getModel", &BaseElementaryMatrix::getModel )
        .def( "getMesh", &BaseElementaryMatrix::getMesh )
        .def( "setMaterialField", &BaseElementaryMatrix::setMaterialField )
        .def( "setElementaryCharacteristics", &BaseElementaryMatrix::setElementaryCharacteristics )
        .def( "setModel", &BaseElementaryMatrix::setModel );

    py::class_< ElementaryMatrixDisplacementReal,
                ElementaryMatrixDisplacementReal::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixDisplacementReal" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementReal, std::string > ) )
        .def( "build", &ElementaryMatrixDisplacementReal::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixDisplacementReal::getFiniteElementDescriptors )
        .def( "getElementaryTerms", &ElementaryMatrixDisplacementReal::getElementaryTerms );

    py::class_< ElementaryMatrixDisplacementComplex,
                ElementaryMatrixDisplacementComplex::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixDisplacementComplex" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementComplex, std::string > ) )
        .def( "build", &ElementaryMatrixDisplacementComplex::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixDisplacementComplex::getFiniteElementDescriptors )
        .def( "getElementaryTerms", &ElementaryMatrixDisplacementComplex::getElementaryTerms );

    py::class_< ElementaryMatrixTemperatureReal,
                ElementaryMatrixTemperatureReal::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixTemperatureReal" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixTemperatureReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixTemperatureReal, std::string > ) )
        .def( "build", &ElementaryMatrixTemperatureReal::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixTemperatureReal::getFiniteElementDescriptors )
        .def( "getElementaryTerms", &ElementaryMatrixTemperatureReal::getElementaryTerms );

    py::class_< ElementaryMatrixPressureComplex,
                ElementaryMatrixPressureComplex::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixPressureComplex" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixPressureComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixPressureComplex, std::string > ) )
        .def( "build", &ElementaryMatrixPressureComplex::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixPressureComplex::getFiniteElementDescriptors )
        .def( "getElementaryTerms", &ElementaryMatrixPressureComplex::getElementaryTerms );
};
