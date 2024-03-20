/**
 * @file DirichletBCInterface.cxx
 * @brief Interface python de DirichletBC
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "aster_pybind.h"

void exportDirichletBCToPython( py::module_ &mod ) {

    py::class_< DirichletBC, DirichletBC::DirichletBCPtr, DataStructure >( mod, "DirichletBC" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def( "build", &DirichletBC::build )
        .def( "getModel", &DirichletBC::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )" )
        .def( "getPhysics", &DirichletBC::getPhysics, R"(
To know the physics supported by the model

Returns:
    str: Mechanics or Thermal or Acoustic
        )" )
        .def( "setSyntax", &DirichletBC::setSyntax, R"(
Function to set the syntax used to build object

Arguments:
    syntax: the syntax
        )",
              py::arg( "syntax" ) );

    py::class_< MechanicalDirichletBC, MechanicalDirichletBC::MechanicalDirichletBCPtr,
                DirichletBC >( mod, "MechanicalDirichletBC" )
        .def( py::init( &initFactoryPtr< MechanicalDirichletBC, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalDirichletBC, std::string, ModelPtr > ) )
        .def( "addBCOnCells",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const std::string & >( &MechanicalDirichletBC::addBCOnCells ) )
        .def( "addBCOnCells",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const VectorString & >( &MechanicalDirichletBC::addBCOnCells ) )
        .def( "addBCOnNodes",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const std::string & >( &MechanicalDirichletBC::addBCOnNodes ) )
        .def( "addBCOnNodes",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const VectorString & >( &MechanicalDirichletBC::addBCOnNodes ) );

    py::class_< ThermalDirichletBC, ThermalDirichletBC::ThermalDirichletBCPtr, DirichletBC >(
        mod, "ThermalDirichletBC" )
        .def( py::init( &initFactoryPtr< ThermalDirichletBC, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ThermalDirichletBC, std::string, ModelPtr > ) )
        .def( "addBCOnCells",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const std::string & >( &ThermalDirichletBC::addBCOnCells ) )
        .def( "addBCOnCells",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const VectorString & >( &ThermalDirichletBC::addBCOnCells ) )
        .def( "addBCOnNodes",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const std::string & >( &ThermalDirichletBC::addBCOnNodes ) )
        .def( "addBCOnNodes",
              py::overload_cast< const PhysicalQuantityComponent &, const ASTERDOUBLE &,
                                 const VectorString & >( &ThermalDirichletBC::addBCOnNodes ) )
        .def( "addBCOnNodes",
              py::overload_cast< const PhysicalQuantityComponent &, const FunctionPtr &,
                                 const VectorString & >( &ThermalDirichletBC::addBCOnNodes ) );

    py::class_< AcousticDirichletBC, AcousticDirichletBC::AcousticDirichletBCPtr, DirichletBC >(
        mod, "AcousticDirichletBC" )
        .def( py::init( &initFactoryPtr< AcousticDirichletBC, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< AcousticDirichletBC, std::string, ModelPtr > ) );
};
