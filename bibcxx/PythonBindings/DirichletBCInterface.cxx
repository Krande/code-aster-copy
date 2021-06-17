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

    bool ( MechanicalDirichletBC::*c1 )( const PhysicalQuantityComponent &,
                                                    const ASTERDOUBLE &, const std::string & ) =
        &MechanicalDirichletBC::addBCOnCells;
    bool ( MechanicalDirichletBC::*c2 )(
        const PhysicalQuantityComponent &, const ASTERDOUBLE &, const VectorString & ) =
        &MechanicalDirichletBC::addBCOnCells;

    bool ( MechanicalDirichletBC::*c3 )( const PhysicalQuantityComponent &,
                                                    const ASTERDOUBLE &, const std::string & ) =
        &MechanicalDirichletBC::addBCOnNodes;
    bool ( MechanicalDirichletBC::*c4 )(
        const PhysicalQuantityComponent &, const ASTERDOUBLE &, const VectorString & ) =
        &MechanicalDirichletBC::addBCOnNodes;

    bool ( ThermalDirichletBC::*c5 )( const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                                 const std::string & ) =
        &ThermalDirichletBC::addBCOnCells;
    bool ( ThermalDirichletBC::*c6 )( const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                                 const VectorString & ) =
        &ThermalDirichletBC::addBCOnCells;

    bool ( ThermalDirichletBC::*c7 )( const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                                 const std::string & ) =
        &ThermalDirichletBC::addBCOnNodes;
    bool ( ThermalDirichletBC::*c8 )( const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                                 const VectorString & ) =
        &ThermalDirichletBC::addBCOnNodes;
    bool ( ThermalDirichletBC::*c9 )( const PhysicalQuantityComponent &,
                                                 const FunctionPtr &,
                                                 const VectorString & ) =
        &ThermalDirichletBC::addBCOnNodes;

    py::class_< DirichletBC, DirichletBC::DirichletBCPtr,
                py::bases< DataStructure > >( "DirichletBC", py::no_init )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "build", &DirichletBC::build )
        .def( "getModel", &DirichletBC::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" ) )  )
        .def( "getPhysics", &DirichletBC::getPhysics, R"(
To know the physics supported by the model

Returns:
    str: Mechanics or Thermal or Acoustic
        )",
              ( py::arg( "self" ) )  );

    py::class_< MechanicalDirichletBC,
                MechanicalDirichletBC::MechanicalDirichletBCPtr,
                py::bases< DirichletBC > >( "MechanicalDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< MechanicalDirichletBC, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< MechanicalDirichletBC, std::string, ModelPtr >))
        .def( "addBCOnCells", c1 )
        .def( "addBCOnCells", c2 )
        .def( "addBCOnNodes", c3 )
        .def( "addBCOnNodes", c4 );

    py::class_< ThermalDirichletBC,
                ThermalDirichletBC::ThermalDirichletBCPtr,
                py::bases< DirichletBC > >( "ThermalDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< ThermalDirichletBC, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< ThermalDirichletBC, std::string, ModelPtr >))
        .def( "addBCOnCells", c5 )
        .def( "addBCOnCells", c6 )
        .def( "addBCOnNodes", c7 )
        .def( "addBCOnNodes", c8 )
        .def( "addBCOnNodes", c9 );

    py::class_< AcousticDirichletBC,
                AcousticDirichletBC::AcousticDirichletBCPtr,
                py::bases< DirichletBC > >( "AcousticDirichletBC", py::no_init )
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< AcousticDirichletBC, ModelPtr >))
        .def( "__init__", py::make_constructor(
            &initFactoryPtr< AcousticDirichletBC, std::string, ModelPtr >));
};
