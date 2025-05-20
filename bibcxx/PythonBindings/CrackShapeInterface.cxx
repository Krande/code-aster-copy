/**
 * @file CrackShapeInterface.cxx
 * @brief Interface python de CrackShape
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "PythonBindings/CrackShapeInterface.h"

#include "aster_pybind.h"

void exportCrackShapeToPython( py::module_ &mod ) {

    py::class_< CrackShape, CrackShape::CrackShapePtr >( mod, "CrackShape" )
        .def( py::init( &initFactoryPtr< CrackShape > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( define_pickling< CrackShape >() )

        .def( "setEllipseCrackShape", &CrackShape::setEllipseCrackShape )
        .def( "setSquareCrackShape", &CrackShape::setSquareCrackShape )
        .def( "setCylinderCrackShape", &CrackShape::setCylinderCrackShape )
        .def( "setNotchCrackShape", &CrackShape::setNotchCrackShape )
        .def( "setHalfPlaneCrackShape", &CrackShape::setHalfPlaneCrackShape )
        .def( "setSegmentCrackShape", &CrackShape::setSegmentCrackShape )
        .def( "setHalfLineCrackShape", &CrackShape::setHalfLineCrackShape )
        .def( "setLineCrackShape", &CrackShape::setLineCrackShape )
        .def( "getShape", &CrackShape::getShape )
        .def( "getShapeName", &CrackShape::getShapeName )
        .def( "getSemiMajorAxis", &CrackShape::getSemiMajorAxis )
        .def( "getSemiMinorAxis", &CrackShape::getSemiMinorAxis )
        .def( "getCenter", &CrackShape::getCenter )
        .def( "getVectX", &CrackShape::getVectX )
        .def( "getVectY", &CrackShape::getVectY )
        .def( "getCrackSide", &CrackShape::getCrackSide )
        .def( "getFilletRadius", &CrackShape::getFilletRadius )
        .def( "getHalfLength", &CrackShape::getHalfLength )
        .def( "getEndPoint", &CrackShape::getEndPoint )
        .def( "getNormal", &CrackShape::getNormal )
        .def( "getTangent", &CrackShape::getTangent )
        .def( "getStartingPoint", &CrackShape::getStartingPoint );
};
