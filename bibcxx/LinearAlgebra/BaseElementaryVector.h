/**
 * @file ElementaryVector.h
 * @brief Definition of elementary vectors
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

#pragma once

#include "astercxx.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCompute.h"
#include "Loads/ListOfLoads.h"
#include "Numbering/DOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Base class for sd_vect_elem
 */
class BaseElementaryVector : public DSWithCppPickling {
  protected:
    /** @brief Flag for empty datastructure (either built or empty)*/
    bool _isBuilt;

    /** @brief Elementary compute */
    ElementaryComputePtr _elemComp;

  public:
    /** @brief Constructor with a name */
    BaseElementaryVector( const std::string name, const std::string type, const ModelPtr model )
        : DSWithCppPickling( name, 19, type ),
          _isBuilt( false ),
          _elemComp( std::make_shared< ElementaryCompute >( getName() ) ) {
        if ( model ) {
            _elemComp->createDescriptor( model );
        }
    };

    /** @brief Constructor with automatic name */
    BaseElementaryVector( const ModelPtr model )
        : BaseElementaryVector( ResultNaming::getNewResultName(), "VECT_ELEM", model ) {};

    py::tuple _getState() const { return py::make_tuple( this->getName(), this->getModel() ); };

  public:
    /**
     * @brief Assembly with dofNume and time (for load)
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                   const ListOfLoadsPtr &loads,
                                                   const ASTERDOUBLE &time = 0. );

    /** @brief Get the model */
    ModelPtr getModel() const { return _elemComp->getModel(); };

    /**
     * @brief Assembly with dofNume and mask on cells
     * @param dofNume object DOFNumbering
     * @param maskCell FieldOnCells to apply mask on each cell
     * @param maskInve flag to inverse mask
     */
    FieldOnNodesRealPtr assembleWithMask( const BaseDOFNumberingPtr &dofNume,
                                          const FieldOnCellsLongPtr &maskCell,
                                          const int &maskInve );

    /**
     * @brief Detect state of datastructure
     * @return true if datastructure has been built (not empty)
     */
    bool isBuilt() { return _isBuilt; };

    /**
     * @brief Set state of datastructure
     * @param bBuilt flag for state of datastructure
     */
    void isBuilt( bool bBuilt ) { _isBuilt = bBuilt; };

    virtual bool build( std::vector< FiniteElementDescriptorPtr > FED = {} ) {
        AS_ABORT( "Not implemented" );
        return false;
    };

    void addSubstructuring( const std::map< std::string, VectorString > &list_load );
};

/**
 * @typedef BaseElementaryVectorPtr
 */
using BaseElementaryVectorPtr = std::shared_ptr< BaseElementaryVector >;
