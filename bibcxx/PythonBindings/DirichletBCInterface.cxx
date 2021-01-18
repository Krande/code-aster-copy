/**
 * @file DirichletBCInterface.cxx
 * @brief Interface python de DirichletBC
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

#include "PythonBindings/DirichletBCInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportDirichletBCToPython() {

    bool ( MechanicalDirichletBCClass::*c1 )( const PhysicalQuantityComponent &,
                                                    const double &, const std::string & ) =
        &MechanicalDirichletBCClass::addBCOnCells;
    bool ( MechanicalDirichletBCClass::*c2 )(
        const PhysicalQuantityComponent &, const double &, const VectorString & ) =
        &MechanicalDirichletBCClass::addBCOnCells;

    bool ( MechanicalDirichletBCClass::*c3 )( const PhysicalQuantityComponent &,
                                                    const double &, const std::string & ) =
        &MechanicalDirichletBCClass::addBCOnNodes;
    bool ( MechanicalDirichletBCClass::*c4 )(
        const PhysicalQuantityComponent &, const double &, const VectorString & ) =
        &MechanicalDirichletBCClass::addBCOnNodes;

    bool ( ThermalDirichletBCClass::*c5 )( const PhysicalQuantityComponent &, const double &,
                                                 const std::string & ) =
        &ThermalDirichletBCClass::addBCOnCells;
    bool ( ThermalDirichletBCClass::*c6 )( const PhysicalQuantityComponent &, const double &,
                                                 const VectorString & ) =
        &ThermalDirichletBCClass::addBCOnCells;

    bool ( ThermalDirichletBCClass::*c7 )( const PhysicalQuantityComponent &, const double &,
                                                 const std::string & ) =
        &ThermalDirichletBCClass::addBCOnNodes;
    bool ( ThermalDirichletBCClass::*c8 )( const PhysicalQuantityComponent &, const double &,
                                                 const VectorString & ) =
        &ThermalDirichletBCClass::addBCOnNodes;
    bool ( ThermalDirichletBCClass::*c9 )( const PhysicalQuantityComponent &,
                                                 const FunctionPtr &,
                                                 const VectorString & ) =
        &ThermalDirichletBCClass::addBCOnNodes;

    py::class_< DirichletBCClass, DirichletBCClass::DirichletBCPtr,
                py::bases< DataStructure > >( "DirichletBC", py::no_init )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "build", &DirichletBCClass::build )
        .def( "getModel", &DirichletBCClass::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" ) )  );

    py::class_< MechanicalDirichletBCClass,
                MechanicalDirichletBCClass::MechanicalDirichletBCPtr,
                py::bases< DirichletBCClass > >( "MechanicalDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< MechanicalDirichletBCClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< MechanicalDirichletBCClass, std::string, ModelPtr >))
        .def( "addBCOnCells", c1 )
        .def( "addBCOnCells", c2 )
        .def( "addBCOnNodes", c3 )
        .def( "addBCOnNodes", c4 );

    py::class_< ThermalDirichletBCClass,
                ThermalDirichletBCClass::ThermalDirichletBCPtr,
                py::bases< DirichletBCClass > >( "ThermalDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< ThermalDirichletBCClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< ThermalDirichletBCClass, std::string, ModelPtr >))
        .def( "addBCOnCells", c5 )
        .def( "addBCOnCells", c6 )
        .def( "addBCOnNodes", c7 )
        .def( "addBCOnNodes", c8 )
        .def( "addBCOnNodes", c9 );

    py::class_< AcousticDirichletBCClass,
                AcousticDirichletBCClass::AcousticDirichletBCPtr,
                py::bases< DirichletBCClass > >( "AcousticDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< AcousticDirichletBCClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< AcousticDirichletBCClass, std::string, ModelPtr >));
};
