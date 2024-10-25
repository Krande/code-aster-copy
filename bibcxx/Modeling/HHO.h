/**
 * @file HHO.h
 * @brief Header of class HHO
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

#pragma once

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Functions/Function.h"
#include "Studies/PhysicalProblem.h"

/**
 * @class HHO
 * @brief Post-processing tools
 */
class HHO {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

    FunctionPtr _createFunc( const ASTERDOUBLE &value ) const;

    /**
     * @brief Project function on HHO space
     */
    FieldOnNodesRealPtr _projectOnHHOSpace( bool faces, const GenericFunctionPtr fct,
                                            ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesRealPtr _projectOnHHOSpace( bool faces, const std::vector< GenericFunctionPtr > fct,
                                            ASTERDOUBLE time = 0.0 ) const;

  public:
    /** @typedef HHOPtr */
    typedef std::shared_ptr< HHO > HHOPtr;

    /** @brief Default constructor disabled */
    HHO( void ) = delete;

    /**
     * @brief Constructor
     * @param PhysicalProblemPtr study
     */
    HHO( const PhysicalProblemPtr &currPhysProblem ) : _phys_problem( currPhysProblem ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    HHO( const py::tuple &tup ) : HHO( tup[0].cast< PhysicalProblemPtr >() ) {};
    py::tuple _getState() const { return py::make_tuple( _phys_problem ); };

    ModelPtr getModel() const;

    /**
     * @brief Project HHO field to H^1-field
     */
    FieldOnNodesRealPtr projectOnLagrangeSpace( const FieldOnNodesRealPtr hho_field ) const;

    /**
     * @brief Evaluate HHO field at quadrature points
     */
    FieldOnCellsRealPtr evaluateAtQuadraturePoints( const FieldOnNodesRealPtr hho_field ) const;

    /**
     * @brief Project H^1-field to HHO field
     */
    FieldOnNodesRealPtr projectOnHHOSpace( const FieldOnNodesRealPtr h1_field ) const;

    /**
     * @brief Project ELGA-field to HHO field
     */
    FieldOnNodesRealPtr projectOnHHOCellSpace( const FieldOnCellsRealPtr field_elga ) const;

    /**
     * @brief Project real value on HHO space
     */
    FieldOnNodesRealPtr projectOnHHOSpace( const ASTERDOUBLE &value ) const;
    FieldOnNodesRealPtr projectOnHHOSpace( const VectorReal &values ) const;

    /**
     * @brief Project real on HHO Cell-space
     */
    FieldOnNodesRealPtr projectOnHHOCellSpace( const ASTERDOUBLE &value ) const;
    FieldOnNodesRealPtr projectOnHHOCellSpace( const VectorReal &values ) const;

    /**
     * @brief Project function on HHO space
     */
    FieldOnNodesRealPtr projectOnHHOSpace( const GenericFunctionPtr fct,
                                           ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesRealPtr projectOnHHOSpace( const std::vector< GenericFunctionPtr > fct,
                                           ASTERDOUBLE time = 0.0 ) const;
    /**
     * @brief Project function on HHO Cell-space
     */
    FieldOnNodesRealPtr projectOnHHOCellSpace( const GenericFunctionPtr fct,
                                               ASTERDOUBLE time = 0.0 ) const;
    FieldOnNodesRealPtr projectOnHHOCellSpace( const std::vector< GenericFunctionPtr > fct,
                                               ASTERDOUBLE time = 0.0 ) const;
};

using HHOPtr = std::shared_ptr< HHO >;
