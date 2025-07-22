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
#include "Discretization/ElementaryCharacteristics.h"
#include "Discretization/ElementaryCompute.h"
#include "Loads/ListOfLoads.h"
#include "Materials/MaterialField.h"
#include "Numbering/DOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Base class for sd_vect_elem
 */
class BaseElementaryVector : public DataStructure {
  protected:
    /** @brief Flag for empty datastructure (either built or empty)*/
    bool _isBuilt;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Elementary compute */
    ElementaryComputePtr _elemComp;

  public:
    /** @brief Constructor with a name */
    BaseElementaryVector( const std::string name, const std::string type = "VECT_ELEM" )
        : DataStructure( name, 19, type ),
          _isBuilt( false ),
          _model( nullptr ),
          _elemComp( std::make_shared< ElementaryCompute >( getName() ) ) {};

    /** @brief Constructor with automatic name */
    BaseElementaryVector() : BaseElementaryVector( ResultNaming::getNewResultName() ) {};

    /** @brief Constructor with automatic name */
    BaseElementaryVector( const ModelPtr model ) : BaseElementaryVector() {
        this->setModel( model );
    };

  public:
    /** @brief Set type of vector */
    void setType( const std::string newType ) { DataStructure::setType( newType ); };

    /**
     * @brief Assembly with dofNume and time (for load)
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                   const ListOfLoadsPtr &loads,
                                                   const ASTERDOUBLE &time = 0. );

    /** @brief Get the model */
    ModelPtr getModel() const { return _model; };

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

    /**
     * @brief Set the model
     * @param currModel pointer to model
     */
    void setModel( const ModelPtr &currModel ) { _model = currModel; };

    /** @brief  Prepare compute */
    void prepareCompute( const std::string option ) { _elemComp->createDescriptor( _model ); };

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
