/**
 * @file ElementaryVectorInterface.cxx
 * @brief Interface python de ElementaryVector
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

#include "PythonBindings/ElementaryVectorInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/FieldOnCellsInterface.h"
#include "PythonBindings/PhysicalProblemInterface.h"

void exportElementaryVectorToPython( py::module_ &mod ) {

    py::class_< BaseElementaryVector, BaseElementaryVectorPtr, DataStructure >(
        mod, "BaseElementaryVector" )
        .def( py::init( &initFactoryPtr< BaseElementaryVector > ) )
        .def( py::init( &initFactoryPtr< BaseElementaryVector, std::string > ) )
        .def( py::init( &initFactoryPtr< BaseElementaryVector, ModelPtr, MaterialFieldPtr,
                                         ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) )
        .def( "addLoad", &BaseElementaryVector::addLoad< const MechanicalLoadRealPtr & > )
        .def( "assembleWithLoadFunctions", &BaseElementaryVector::assembleWithLoadFunctions,
              py::arg( "dofNume" ), py::arg( "time" ) = 0. )
        .def( "assembleWithMask", &BaseElementaryVector::assembleWithMask )
        .def( "setType", &BaseElementaryVector::setType )
        .def( "setListOfLoads", &BaseElementaryVector::setListOfLoads )
        .def( "setMaterialField", &BaseElementaryVector::setMaterialField )
        .def( "setModel", &BaseElementaryVector::setModel )
        .def( "prepareCompute", &BaseElementaryVector::prepareCompute )
        .def( "setElementaryCharacteristics", &BaseElementaryVector::setElementaryCharacteristics )
        .def( "addSubstructuring", &BaseElementaryVector::addSubstructuring )
        .def( "build", &BaseElementaryVector::build,
              py::arg( "FED" ) = std::vector< FiniteElementDescriptorPtr >() );

    py::class_< ElementaryVectorReal, ElementaryVectorRealPtr, BaseElementaryVector >(
        mod, "ElementaryVectorReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal, ModelPtr, MaterialFieldPtr,
                                         ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) )
        .def( "getVeass", &ElementaryVectorReal::getVeass )
        .def( "assemble", &ElementaryVectorReal::assemble )
        .def( "addElementaryTerm",
              py::overload_cast< const ElementaryTermRealPtr & >(
                  &ElementaryVectorReal::addElementaryTerm ),
              R"(
            Add elementary term

            Arguments:
                term (ElementaryTermReal): elementary term
            )",
              py::arg( "term" ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermRealPtr > & >(
                  &ElementaryVectorReal::addElementaryTerm ),
              R"(
            Add vector of elementary term

            Arguments:
                terms (list[ElementaryTermReal]): vector of elementary term
            )",
              py::arg( "terms" ) )
        .def( "getElementaryTerms", &ElementaryVectorReal::getElementaryTerms )
        .def( "getFiniteElementDescriptor", &ElementaryVectorReal::getFiniteElementDescriptor );

    py::class_< ElementaryVectorComplex, ElementaryVectorComplexPtr, BaseElementaryVector >(
        mod, "ElementaryVectorComplex" )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex, ModelPtr, MaterialFieldPtr,
                                         ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) )
        .def( "getVeass", &ElementaryVectorComplex::getVeass )
        .def( "addElementaryTerm",
              py::overload_cast< const ElementaryTermComplexPtr & >(
                  &ElementaryVectorComplex::addElementaryTerm ),
              R"(
            Add elementary term

            Arguments:
                term (ElementaryTermComplex): elementary term
            )",
              py::arg( "term" ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermComplexPtr > & >(
                  &ElementaryVectorComplex::addElementaryTerm ),
              R"(
            Add vector of elementary term

            Arguments:
                terms (list[ElementaryTermComplex]): vector of elementary term
            )",
              py::arg( "terms" ) )
        .def( "getElementaryTerms", &ElementaryVectorComplex::getElementaryTerms )
        .def( "assemble", &ElementaryVectorComplex::assemble );

    py::class_< ElementaryVectorDisplacementReal, ElementaryVectorDisplacementRealPtr,
                ElementaryVectorReal >( mod, "ElementaryVectorDisplacementReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorDisplacementReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorDisplacementReal, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryVectorDisplacementReal, ModelPtr, MaterialFieldPtr,
                                       ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) );

    py::class_< ElementaryVectorTemperatureReal, ElementaryVectorTemperatureRealPtr,
                ElementaryVectorReal >( mod, "ElementaryVectorTemperatureReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorTemperatureReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorTemperatureReal, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryVectorTemperatureReal, ModelPtr, MaterialFieldPtr,
                                       ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) );

    py::class_< ElementaryVectorPressureComplex, ElementaryVectorPressureComplexPtr,
                ElementaryVectorComplex >( mod, "ElementaryVectorPressureComplex" )
        .def( py::init( &initFactoryPtr< ElementaryVectorPressureComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorPressureComplex, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryVectorPressureComplex, ModelPtr, MaterialFieldPtr,
                                       ElementaryCharacteristicsPtr, ListOfLoadsPtr > ) );
};
