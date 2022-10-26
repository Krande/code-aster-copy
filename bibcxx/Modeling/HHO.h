
/**
 * @file HHO.h
 * @brief Header of class HHO
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

#include "astercxx.h"

#include "Behaviours/BehaviourProperty.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Numbering/DOFNumbering.h"
#include "Studies/PhysicalProblem.h"

/**
 * @class HHO
 * @brief Post-processing tools
 */
class HHO {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

  public:
    /** @typedef HHOPtr */
    typedef std::shared_ptr< HHO > HHOPtr;

    /** @brief Default constructor disabled */
    HHO( void ) = delete;

    /**
     * @brief Constructor
     * @param PhysicalProblemPtr study
     */
    HHO( const PhysicalProblemPtr &currPhysProblem ) : _phys_problem( currPhysProblem ){};

    /**
     * @brief Project HHO field to H^1-field
     */
    FieldOnNodesRealPtr projectOnLagrangeSpace( const FieldOnNodesRealPtr hho_field ) const;

    /**
     * @brief Project H^1-field on HHO space
     */
    FieldOnNodesRealPtr projectOnHHOSpace( const ASTERDOUBLE &value ) const;
};

using HHOPtr = std::shared_ptr< HHO >;
