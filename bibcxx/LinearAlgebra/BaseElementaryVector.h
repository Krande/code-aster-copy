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

#include "astercxx.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Discretization/ElementaryCompute.h"
#include "Loads/ListOfLoads.h"
#include "Materials/MaterialField.h"
#include "Numbering/DOFNumbering.h"
#include "Studies/PhysicalProblem.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Base class for sd_vect_elem
 */
class BaseElementaryVector : public DataStructure {
  protected:
    /** @brief Flag for empty datastructure */
    bool _isEmpty;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Field of material parameters */
    MaterialFieldPtr _materialField;

    /** @brief Elementary characteristics */
    ElementaryCharacteristicsPtr _elemChara;

    /** @brief Liste de charges */
    ListOfLoadsPtr _listOfLoads;

    /** @brief Elementary compute */
    ElementaryComputePtr _elemComp;

    /**
     * @brief Set the option
     * @param currOption option
     */
    void setOption( const std::string option ) { _elemComp->setOption( option ); };

  public:
    /** @brief Constructor with a name */
    BaseElementaryVector( const std::string name, const std::string type = "VECT_ELEM" )
        : DataStructure( name, 19, type ),
          _isEmpty( true ),
          _model( nullptr ),
          _materialField( nullptr ),
          _elemChara( nullptr ),
          _elemComp( std::make_shared< ElementaryCompute >( getName() ) ),
          _listOfLoads( std::make_shared< ListOfLoads >() ){};

    /** @brief Constructor with automatic name */
    BaseElementaryVector() : BaseElementaryVector( ResultNaming::getNewResultName() ){};

    /** @brief Constructor with automatic name */
    BaseElementaryVector( const PhysicalProblemPtr phys_pb ) : BaseElementaryVector() {
        this->setPhysicalProblem( phys_pb );
    };

  public:
    /** @brief Add a load to elementary vector */
    template < typename... Args >
    void addLoad( const Args &... a ) {
        _listOfLoads->addLoad( a... );
    };

    /** @brief Set type of vector */
    void setType( const std::string newType ) { DataStructure::setType( newType ); };

    /**
     * @brief Assembly with dofNume and time (for load)
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                   const ASTERDOUBLE &time = 0. );

    /** @brief Get the model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get option */
    std::string getOption() const { return _elemComp->getOption(); };

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
     * @return true if empty datastructure
     */
    bool isEmpty() { return _isEmpty; };

    /**
     * @brief Set state of datastructure
     * @param bEmpty flag for state of datastructure
     */
    void isEmpty( bool bEmpty ) { _isEmpty = bEmpty; };

    /**
     * @brief Set physical problem
     */
    void setPhysicalProblem( const PhysicalProblemPtr phys_pb ) {
        _model = phys_pb->getModel();
        _materialField = phys_pb->getMaterialField();
        _elemChara = phys_pb->getElementaryCharacteristics();
        _listOfLoads = phys_pb->getListOfLoads();
    };

    /**
     * @brief Set list of loads
     * @param currentList list of loads
     */
    void setListOfLoads( const ListOfLoadsPtr &currentList ) { _listOfLoads = currentList; };

    /**
     * @brief Set the field of material parameters
     * @param currMaterialField pointer to material field
     */
    void setMaterialField( const MaterialFieldPtr &currMaterialField ) {
        _materialField = currMaterialField;
    };

    /**
     * @brief Set elementary characteristics
     * @param currElemChara pointer to elementary characteristics
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &currElemChara ) {
        _elemChara = currElemChara;
    };

    /**
     * @brief Set the model
     * @param currModel pointer to model
     */
    void setModel( const ModelPtr &currModel ) { _model = currModel; };

    /** @brief  Prepare compute */
    void prepareCompute( const std::string option ) {
        setOption( option );
        _elemComp->setOption( option );
        if ( option != "WRAP_FORTRAN" ) {
            _elemComp->createDescriptor( _model, _materialField, _elemChara );
        }
    };
};

/**
 * @typedef BaseElementaryVectorPtr
 */
using BaseElementaryVectorPtr = std::shared_ptr< BaseElementaryVector >;
