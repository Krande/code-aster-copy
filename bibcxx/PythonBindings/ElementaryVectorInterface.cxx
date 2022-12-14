/**
 * @file ElementaryVectorInterface.cxx
 * @brief Interface python de ElementaryVector
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

#include "PythonBindings/ElementaryVectorInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/FieldOnCellsInterface.h"
#include "PythonBindings/PhysicalProblemInterface.h"

void exportElementaryVectorToPython( py::module_ &mod ) {

    py::class_< BaseElementaryVector, BaseElementaryVectorPtr, DataStructure >(
        mod, "BaseElementaryVector" )
        .def( py::init( &initFactoryPtr< BaseElementaryVector > ) )
        .def( py::init( &initFactoryPtr< BaseElementaryVector, std::string > ) )
        .def( py::init( &initFactoryPtr< BaseElementaryVector, PhysicalProblemPtr > ) )
        .def( "addLoad", &BaseElementaryVector::addLoad< const MechanicalLoadRealPtr & > )
        .def( "assembleWithLoadFunctions", &BaseElementaryVector::assembleWithLoadFunctions,
              py::arg( "dofNume" ), py::arg( "time" ) = 0. )
        .def( "assembleWithMask", &BaseElementaryVector::assembleWithMask )
        .def( "setType", &BaseElementaryVector::setType )
        .def( "setPhysicalProblem", &BaseElementaryVector::setPhysicalProblem, R"(
            Set the physical problem

            Arguments:
                phys_pb (PhysicalProblem): the physical problem.
            )",
              py::arg( "phys_pb" ) )
        .def( "setListOfLoads", &BaseElementaryVector::setListOfLoads )
        .def( "setMaterialField", &BaseElementaryVector::setMaterialField )
        .def( "setModel", &BaseElementaryVector::setModel )
        .def( "setElementaryCharacteristics", &BaseElementaryVector::setElementaryCharacteristics )
        .def( "build", &BaseElementaryVector::build );

    py::class_< ElementaryVectorReal, ElementaryVectorRealPtr, BaseElementaryVector >(
        mod, "ElementaryVectorReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorReal, PhysicalProblemPtr > ) )
        .def( "getVeass", &ElementaryVectorReal::getVeass )
        .def( "assemble", &ElementaryVectorReal::assemble );

    py::class_< ElementaryVectorComplex, ElementaryVectorComplexPtr, BaseElementaryVector >(
        mod, "ElementaryVectorComplex" )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorComplex, PhysicalProblemPtr > ) )
        .def( "getVeass", &ElementaryVectorComplex::getVeass )
        .def( "assemble", &ElementaryVectorComplex::assemble );

    py::class_< ElementaryVectorDisplacementReal, ElementaryVectorDisplacementRealPtr,
                ElementaryVectorReal >( mod, "ElementaryVectorDisplacementReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorDisplacementReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorDisplacementReal, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorDisplacementReal, PhysicalProblemPtr > ) );

    py::class_< ElementaryVectorTemperatureReal, ElementaryVectorTemperatureRealPtr,
                ElementaryVectorReal >( mod, "ElementaryVectorTemperatureReal" )
        .def( py::init( &initFactoryPtr< ElementaryVectorTemperatureReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorTemperatureReal, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorTemperatureReal, PhysicalProblemPtr > ) );

    py::class_< ElementaryVectorPressureComplex, ElementaryVectorPressureComplexPtr,
                ElementaryVectorComplex >( mod, "ElementaryVectorPressureComplex" )
        .def( py::init( &initFactoryPtr< ElementaryVectorPressureComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorPressureComplex, std::string > ) )
        .def( py::init( &initFactoryPtr< ElementaryVectorPressureComplex, PhysicalProblemPtr > ) );
};
