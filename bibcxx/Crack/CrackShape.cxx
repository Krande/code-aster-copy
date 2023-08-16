/**
 * @file CrackShape.cxx
 * @brief Class to describe the possible shape of cracks for XFEM
 * @author Nicolas Tardieu
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
#include "CrackShape.h"

CrackShape::CrackShape() { _shape = Shape::NoShape; };

py::tuple CrackShape::_getState() const {
    return py::make_tuple( _semiMajorAxis, _semiMinorAxis, _center, _vectX, _vectY, _crackSide,
                           _filletRadius, _halfLength, _endPoint, _normal, _tangent,
                           _startingPoint );
}

CrackShape::CrackShape( const py::tuple &tup ) : CrackShape() {
    _shape = tup[0].cast< Shape >();
    _semiMajorAxis = tup[1].cast< ASTERDOUBLE >();
    _semiMinorAxis = tup[2].cast< ASTERDOUBLE >();
    _center = tup[3].cast< VectorReal >();
    _vectX = tup[4].cast< VectorReal >();
    _vectY = tup[5].cast< VectorReal >();
    _crackSide = tup[6].cast< std::string >();
    _filletRadius = tup[7].cast< ASTERDOUBLE >();
    _halfLength = tup[8].cast< ASTERDOUBLE >();
    _endPoint = tup[9].cast< VectorReal >();
    _normal = tup[10].cast< VectorReal >();
    _tangent = tup[11].cast< VectorReal >();
    _startingPoint = tup[12].cast< VectorReal >();
}

void CrackShape::setEllipseCrackShape( ASTERDOUBLE semiMajorAxis, ASTERDOUBLE semiMinorAxis,
                                       VectorReal center, VectorReal vectX, VectorReal vectY,
                                       std::string crackSide ) {
    _shape = Shape::Ellipse;
    _semiMajorAxis = semiMajorAxis;
    _semiMinorAxis = semiMinorAxis;
    _center = center;
    _vectX = vectX;
    _vectY = vectY;
    _crackSide = crackSide;
};

void CrackShape::setSquareCrackShape( ASTERDOUBLE semiMajorAxis, ASTERDOUBLE semiMinorAxis,
                                      ASTERDOUBLE filletRadius, VectorReal center, VectorReal vectX,
                                      VectorReal vectY, std::string crackSide ) {
    _shape = Shape::Square;
    _semiMajorAxis = semiMajorAxis;
    _semiMinorAxis = semiMinorAxis;
    _filletRadius = filletRadius;
    _center = center;
    _vectX = vectX;
    _vectY = vectY;
    _crackSide = crackSide;
};

void CrackShape::setCylinderCrackShape( ASTERDOUBLE semiMajorAxis, ASTERDOUBLE semiMinorAxis,
                                        VectorReal center, VectorReal vectX, VectorReal vectY ) {
    _shape = Shape::Cylinder;
    _semiMajorAxis = semiMajorAxis;
    _semiMinorAxis = semiMinorAxis;
    _center = center;
    _vectX = vectX;
    _vectY = vectY;
};

void CrackShape::setNotchCrackShape( ASTERDOUBLE halfLength, ASTERDOUBLE filletRadius,
                                     VectorReal center, VectorReal vectX, VectorReal vectY ) {
    _shape = Shape::Notch;
    _halfLength = halfLength;
    _filletRadius = filletRadius;
    _center = center;
    _vectX = vectX;
    _vectY = vectY;
};

void CrackShape::setHalfPlaneCrackShape( VectorReal endPoint, VectorReal normal,
                                         VectorReal tangent ) {
    _shape = Shape::HalfPlane;
    _endPoint = endPoint;
    _normal = normal;
    _tangent = tangent;
};
void CrackShape::setSegmentCrackShape( VectorReal startingPoint, VectorReal endPoint ) {
    _shape = Shape::Segment;
    _startingPoint = startingPoint;
    _endPoint = endPoint;
};

void CrackShape::setHalfLineCrackShape( VectorReal startingPoint, VectorReal tangent ) {
    _shape = Shape::HalfLine;
    _startingPoint = startingPoint;
    _tangent = tangent;
};

void CrackShape::setLineCrackShape( VectorReal startingPoint, VectorReal tangent ) {
    _shape = Shape::Line;
    _startingPoint = startingPoint;
    _tangent = tangent;
};

std::string CrackShape::getShapeName() const {
    if ( _shape == Shape::NoShape )
        return "NoShape";
    else if ( _shape == Shape::Ellipse )
        return "ELLIPSE";
    else if ( _shape == Shape::Square )
        return "RECTANGLE";
    else if ( _shape == Shape::Cylinder )
        return "CYLINDRE";
    else if ( _shape == Shape::Notch )
        return "ENTAILLE";
    else if ( _shape == Shape::HalfPlane )
        return "DEMI_PLAN";
    else if ( _shape == Shape::Segment )
        return "SEGMENT";
    else if ( _shape == Shape::HalfLine )
        return "DEMI_DROITE";
    else if ( _shape == Shape::Line )
        return "DROITE";
    else {
        throw std::runtime_error( "Unknown shape" );
    }

    return std::string();
};
