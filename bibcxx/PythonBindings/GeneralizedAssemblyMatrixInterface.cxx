/**
 * @file GeneralizedAssemblyMatrixInterface.cxx
 * @brief Interface python de GeneralizedAssemblyMatrix
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

#include "PythonBindings/GeneralizedAssemblyMatrixInterface.h"

#include "aster_pybind.h"

#include "Numbering/GeneralizedDOFNumbering.h"
#include "PythonBindings/VariantModalBasisInterface.h"

void exportGeneralizedAssemblyMatrixToPython( py::module_ &mod ) {

    py::class_< GenericGeneralizedAssemblyMatrix, GenericGeneralizedAssemblyMatrixPtr,
                DataStructure >( mod, "GeneralizedAssemblyMatrix" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "getGeneralizedDOFNumbering",
              &GenericGeneralizedAssemblyMatrix::getGeneralizedDOFNumbering )
        .def( "getModalBasis", &getModalBasis< GenericGeneralizedAssemblyMatrixPtr > )
        .def( "setGeneralizedDOFNumbering",
              &GenericGeneralizedAssemblyMatrix::setGeneralizedDOFNumbering )
        .def( "setModalBasis", py::overload_cast< const GeneralizedModeResultPtr & >(
                                   &GenericGeneralizedAssemblyMatrix::setModalBasis ) )
        .def( "setModalBasis", py::overload_cast< const ModeResultPtr & >(
                                   &GenericGeneralizedAssemblyMatrix::setModalBasis ) );

    py::class_< GeneralizedAssemblyMatrixReal, GeneralizedAssemblyMatrixRealPtr,
                GenericGeneralizedAssemblyMatrix >( mod, "GeneralizedAssemblyMatrixReal" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixReal > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixReal, std::string > ) );

    py::class_< GeneralizedAssemblyMatrixComplex, GeneralizedAssemblyMatrixComplexPtr,
                GenericGeneralizedAssemblyMatrix >( mod, "GeneralizedAssemblyMatrixComplex" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixComplex > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyMatrixComplex, std::string > ) );
};
