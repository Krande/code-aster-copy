/**
 * @file GeneralizedResultInterface.cxx
 * @brief Interface python de GeneralizedResult
 * @author Natacha Béreux
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


#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/GeneralizedResultInterface.h"


void exportGeneralizedResultToPython()
{


    py::class_< GeneralizedResultReal,
            GeneralizedResultRealPtr,
            py::bases< DataStructure > >
            ( "GeneralizedResultReal", py::no_init )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
    ;

    py::class_< GeneralizedResultComplex,
            GeneralizedResultComplexPtr,
            py::bases< DataStructure > >
            ( "GeneralizedResultComplex", py::no_init )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
    ;

    py::class_< TransientGeneralizedResult,
            TransientGeneralizedResultPtr,
            py::bases< GeneralizedResultReal > >
            ( "TransientGeneralizedResult", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< TransientGeneralizedResult >) )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< TransientGeneralizedResult,
                             std::string >) )
        .def( "setGeneralizedDOFNumbering",
              &TransientGeneralizedResult::setGeneralizedDOFNumbering )
        .def( "getGeneralizedDOFNumbering",
              &TransientGeneralizedResult::getGeneralizedDOFNumbering )
        .def( "setDOFNumbering",
              &TransientGeneralizedResult::setDOFNumbering )
        .def( "getDOFNumbering",
              &TransientGeneralizedResult::getDOFNumbering )
    ;

    py::class_< HarmoGeneralizedResult,
            HarmoGeneralizedResultPtr,
            py::bases< GeneralizedResultComplex > >
            ( "HarmoGeneralizedResult", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< HarmoGeneralizedResult >) )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< HarmoGeneralizedResult,
                             std::string >) )
        .def( "getGeneralizedDOFNumbering",
              &HarmoGeneralizedResult::getGeneralizedDOFNumbering )
        .def( "setGeneralizedDOFNumbering",
              &HarmoGeneralizedResult::setGeneralizedDOFNumbering )
        .def( "setDOFNumbering",
              &HarmoGeneralizedResult::setDOFNumbering )
        .def( "getDOFNumbering",
              &HarmoGeneralizedResult::getDOFNumbering )
    ;
};
