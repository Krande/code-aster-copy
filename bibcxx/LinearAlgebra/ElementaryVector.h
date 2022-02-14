#ifndef ELEMENTARYVECTOR_H_
#define ELEMENTARYVECTOR_H_

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

#include "astercxx.h"

#include "DataFields/ElementaryTerm.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Discretization/ElementaryCompute.h"
#include "Loads/ListOfLoads.h"
#include "Loads/PhysicalQuantity.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/DOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Generic class for sd_vect_elem
 */
class BaseElementaryVector : public DataStructure {
  protected:
    /** @brief Option to compute */
    std::string _option;

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
    void setOption( const std::string option ) { _option = option; };

  public:
    /** @brief Constructor with a name */
    BaseElementaryVector( const std::string name, const std::string type = "VECT_ELEM" )
        : DataStructure( name, 19, type ),
          _isEmpty( true ),
          _model( nullptr ),
          _materialField( nullptr ),
          _elemChara( nullptr ),
          _elemComp( nullptr ),
          _option( " " ),
          _listOfLoads( new ListOfLoads() ){};

    /** @brief Constructor with automatic name */
    BaseElementaryVector() : BaseElementaryVector( ResultNaming::getNewResultName() ){};

  public:
    /** @brief Add a load to elementary vector */
    template < typename... Args >
    void addLoad( const Args &...a ) {
        _listOfLoads->addLoad( a... );
    };

    /** @brief Set type of vector */
    void setType( const std::string newType ) { DataStructure::setType( newType ); };

    /**
     * @brief Assembly with dofNume and time (for load)
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume,
                                                   const ASTERDOUBLE &time );

    /**
     * @brief Assembly with dofNume and apply functions from loads
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assembleWithLoadFunctions( const BaseDOFNumberingPtr &dofNume ) {
        return assembleWithLoadFunctions( dofNume, 0. );
    };

    /**
     * @brief Assembly with dofNume
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesRealPtr assemble( const BaseDOFNumberingPtr dofNume ) const;

    /** @brief Get the model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get the mesh */
    BaseMeshPtr getMesh( void ) const;

    /** @brief Get option */
    std::string getOption() const { return _option; };

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
        _elemComp = boost::make_shared< ElementaryCompute >( this->getName(), _option );
        if ( _option != "WRAP_FORTRAN" ) {
            _elemComp->createDescriptor( _model, _materialField, _elemChara );
            _elemComp->createListOfElementaryTerms();
        }
    };

    std::string generateNameOfElementaryTerm() const {
        std::string elemTermName( " " );
        std::ostringstream numString;
        numString << std::setw( 7 ) << std::setfill( '0' ) << _elemComp->getIndexName();
        _elemComp->nextIndexName();
        elemTermName = this->getName().substr( 0, 8 ) + "." + numString.str();
        return elemTermName;
    }
};

/** @typedef BaseElementaryVectorPtr */
typedef boost::shared_ptr< BaseElementaryVector > BaseElementaryVectorPtr;

/**
 * @class ElementaryVector
 * @brief Class for sd_vect_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryVector : public BaseElementaryVector {
  private:
    /** @brief Vectors of RESUELEM */
    std::vector< boost::shared_ptr< ElementaryTerm< ValueType > > > _elemTerm;

  public:
    /** @typedef ElementaryVectorPtr */
    typedef boost::shared_ptr< ElementaryVector< ValueType, PhysicalQuantity > >
        ElementaryVectorPtr;

    /** @brief Constructor with a name */
    ElementaryVector( const std::string name )
        : BaseElementaryVector(
              name, "VECT_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ) {
        _elemTerm.clear();
    };

    /** @brief Constructor with automatic name */
    ElementaryVector() : ElementaryVector( ResultNaming::getNewResultName() ){};

    /**
     * @brief Function to update ElementaryTerm
     */
    bool build() {
        _elemComp = boost::make_shared< ElementaryCompute >( this->getName() );
        if ( _elemComp->hasElementaryTerm() ) {
            std::vector< JeveuxChar24 > elemTermNames = _elemComp->getNameOfElementaryTerms();
            for ( int pos = 0; pos < elemTermNames.size(); ++pos ) {
                const std::string name = elemTermNames[pos].toString();
                if ( trim( name ) != "" ) {
                    boost::shared_ptr< ElementaryTerm< ValueType > > toPush(
                        new ElementaryTerm< ValueType >( name ) );
                    _elemTerm.push_back( toPush );
                }
            }
        }
        return true;
    };

    /**
     * @brief Add elementary term
     */
    void addElementaryTerm( const boost::shared_ptr< ElementaryTerm< ValueType > > &elemTerm ) {
        _elemComp->addElementaryTerm( elemTerm->getName() );
        _elemTerm.push_back( elemTerm );
    };

    friend class DiscreteComputation;
};

/** @typedef Elementary vector for displacement-double */
template class ElementaryVector< ASTERDOUBLE, Displacement >;
typedef ElementaryVector< ASTERDOUBLE, Displacement > ElementaryVectorDisplacementReal;
typedef boost::shared_ptr< ElementaryVectorDisplacementReal > ElementaryVectorDisplacementRealPtr;

/** @typedef Elementary vector for temperature-double */
template class ElementaryVector< ASTERDOUBLE, Temperature >;
typedef ElementaryVector< ASTERDOUBLE, Temperature > ElementaryVectorTemperatureReal;
typedef boost::shared_ptr< ElementaryVectorTemperatureReal > ElementaryVectorTemperatureRealPtr;

/** @typedef Elementary vector for pressure-complex */
template class ElementaryVector< ASTERCOMPLEX, Pressure >;
typedef ElementaryVector< ASTERCOMPLEX, Pressure > ElementaryVectorPressureComplex;
typedef boost::shared_ptr< ElementaryVectorPressureComplex > ElementaryVectorPressureComplexPtr;

#endif /* ELEMENTARYVECTOR_H_ */
