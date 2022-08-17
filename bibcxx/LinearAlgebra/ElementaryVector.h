/**
 * @file ElementaryVector.h
 * @brief Definition of elementary vectors
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

#pragma once

#include "LinearAlgebra/GenericElementaryVector.h"

/**
 * @class ElementaryVector
 * @brief Class for sd_vect_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryVector : public GenericElementaryVector< ValueType > {

  public:
    /** @typedef ElementaryVectorPtr */
    typedef std::shared_ptr< ElementaryVector< ValueType, PhysicalQuantity > > ElementaryVectorPtr;

    /** @brief Constructor with a name */
    ElementaryVector( const std::string name )
        : GenericElementaryVector< ValueType >(
              name, "VECT_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ){};

    /** @brief Constructor with automatic name */
    ElementaryVector() : ElementaryVector( ResultNaming::getNewResultName() ){};

    ElementaryVector( const PhysicalProblemPtr phys_pb ) : ElementaryVector() {
        this->setPhysicalProblem( phys_pb );
    };
};

/** @typedef Elementary vector for displacement-double */
template class ElementaryVector< ASTERDOUBLE, Displacement >;
typedef ElementaryVector< ASTERDOUBLE, Displacement > ElementaryVectorDisplacementReal;
typedef std::shared_ptr< ElementaryVectorDisplacementReal > ElementaryVectorDisplacementRealPtr;

/** @typedef Elementary vector for temperature-double */
template class ElementaryVector< ASTERDOUBLE, Temperature >;
typedef ElementaryVector< ASTERDOUBLE, Temperature > ElementaryVectorTemperatureReal;
typedef std::shared_ptr< ElementaryVectorTemperatureReal > ElementaryVectorTemperatureRealPtr;

/** @typedef Elementary vector for pressure-complex */
template class ElementaryVector< ASTERCOMPLEX, Pressure >;
typedef ElementaryVector< ASTERCOMPLEX, Pressure > ElementaryVectorPressureComplex;
typedef std::shared_ptr< ElementaryVectorPressureComplex > ElementaryVectorPressureComplexPtr;
