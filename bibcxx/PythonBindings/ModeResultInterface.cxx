/**
 * @file ModeResultInterface.cxx
 * @brief Interface python de ModeResult
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

#include "PythonBindings/ModeResultInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/VariantStiffnessMatrixInterface.h"

void exportModeResultToPython() {

    bool ( ModeResult::*c1 )( const AssemblyMatrixDisplacementRealPtr & ) =
        &ModeResult::setStiffnessMatrix;
    bool ( ModeResult::*c2 )( const AssemblyMatrixTemperatureRealPtr & ) =
        &ModeResult::setStiffnessMatrix;
    bool ( ModeResult::*c3 )( const AssemblyMatrixDisplacementComplexPtr & ) =
        &ModeResult::setStiffnessMatrix;
    bool ( ModeResult::*c4 )( const AssemblyMatrixPressureRealPtr & ) =
        &ModeResult::setStiffnessMatrix;

    bool ( ModeResult::*c5 )( const AssemblyMatrixDisplacementRealPtr & ) =
        &ModeResult::setMassMatrix;
    bool ( ModeResult::*c6 )( const AssemblyMatrixTemperatureRealPtr & ) =
        &ModeResult::setMassMatrix;
    bool ( ModeResult::*c7 )( const AssemblyMatrixDisplacementComplexPtr & ) =
        &ModeResult::setMassMatrix;
    bool ( ModeResult::*c8 )( const AssemblyMatrixPressureRealPtr & ) =
        &ModeResult::setMassMatrix;

    py::class_< ModeResult, ModeResultPtr,
                py::bases< FullResult > >( "ModeResult",
                                                             py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< ModeResult >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ModeResult, std::string >))
        .def( "getDOFNumbering", &ModeResult::getDOFNumbering )
        .def("getStiffnessMatrix", &getStiffnessMatrix< ModeResultPtr >)
        .def( "setStiffnessMatrix", c1 )
        .def( "setStiffnessMatrix", c2 )
        .def( "setStiffnessMatrix", c3 )
        .def( "setStiffnessMatrix", c4 )
        .def("getMassMatrix", &getStiffnessMatrix< ModeResultPtr >)
        .def( "setMassMatrix", c5 )
        .def( "setMassMatrix", c6 )
        .def( "setMassMatrix", c7 )
        .def( "setMassMatrix", c8 )
        .def( "setStructureInterface", &ModeResult::setStructureInterface );
};

void exportModeResultComplexToPython() {

    bool ( ModeResultComplex::*c1 )(
        const AssemblyMatrixDisplacementRealPtr & ) =
        &ModeResultComplex::setStiffnessMatrix;
    bool ( ModeResultComplex::*c2 )(
        const AssemblyMatrixDisplacementComplexPtr & ) =
        &ModeResultComplex::setStiffnessMatrix;
    bool ( ModeResultComplex::*c3 )(
        const AssemblyMatrixDisplacementComplexPtr & ) =
        &ModeResultComplex::setStiffnessMatrix;
    bool ( ModeResultComplex::*c4 )(
        const AssemblyMatrixPressureRealPtr & ) =
        &ModeResultComplex::setStiffnessMatrix;

    py::class_< ModeResultComplex, ModeResultComplexPtr,
                py::bases< ModeResult > >( "ModeResultComplex",
                                                                py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ModeResultComplex >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ModeResultComplex, std::string >))
        .def( "setDampingMatrix", &ModeResultComplex::setDampingMatrix )
        .def( "setStiffnessMatrix", c1 )
        .def( "setStiffnessMatrix", c2 )
        .def( "setStiffnessMatrix", c3 )
        .def( "setStiffnessMatrix", c4 )
        .def( "setStructureInterface",
              &ModeResultComplex::setStructureInterface );
};
