#ifndef ELEMENTARYMATRIX_H_
#define ELEMENTARYMATRIX_H_

/**
 * @file ElementaryMatrix.h
 * @brief Definition of elementary matrices
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
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/PhysicalQuantity.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryMatrix
 * @brief Generic class for sd_matr_elem
 */
class BaseElementaryMatrix : public DataStructure {
  protected:
    /** @brief Objects for computation of elementary terms */
    JeveuxVectorChar24 _description;
    JeveuxVectorChar24 _listOfElementaryTerms;

    /** @brief Flag for empty datastructure */
    bool _isEmpty;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Field of material parameters */
    MaterialFieldPtr _materialField;

    /** @brief Elementary characteristics */
    ElementaryCharacteristicsPtr _elemChara;

    /** @brief Vectors of FiniteElementDescriptor */
    std::vector< FiniteElementDescriptorPtr > _FEDVector;
    std::set< std::string > _FEDNames;

  public:
    /** @brief Constructor with a name */
    BaseElementaryMatrix( const std::string name, const std::string type = "MATR_ELEM" )
        : DataStructure( name, 19, type ),
          _description( JeveuxVectorChar24( getName() + ".RERR" ) ),
          _listOfElementaryTerms( JeveuxVectorChar24( getName() + ".RELR" ) ),
          _isEmpty( true ),
          _model( nullptr ),
          _materialField( nullptr ),
          _elemChara( nullptr ){};

    /** @brief Constructor with automatic name */
    BaseElementaryMatrix( const std::string type = "MATR_ELEM" )
        : BaseElementaryMatrix( ResultNaming::getNewResultName(), type ){};

  public:
    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptor( const FiniteElementDescriptorPtr &modelFED ) {
        const auto name = trim( modelFED->getName() );
        if ( _FEDNames.find( name ) == _FEDNames.end() ) {
            _FEDVector.push_back( _model->getFiniteElementDescriptor() );
            _FEDNames.insert( name );
            return true;
        }
        return false;
    };

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() { return _FEDVector; };

    /** @brief Get the field of material parameters */
    MaterialFieldPtr getMaterialField() const {
        if ( _materialField == nullptr )
            throw std::runtime_error( "MaterialField is not set" );
        return _materialField;
    };

    /** @brief Get the model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get the mesh */
    BaseMeshPtr getMesh( void ) const;

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
    void setModel( const ModelPtr &currModel ) {
        _model = currModel;
        auto modelFED = _model->getFiniteElementDescriptor();
        const auto name = trim( modelFED->getName() );
        if ( _FEDNames.find( name ) == _FEDNames.end() ) {
            _FEDVector.push_back( modelFED );
            _FEDNames.insert( name );
        }
    };
};

/** @typedef BaseElementaryMatrixPtr */
typedef boost::shared_ptr< BaseElementaryMatrix > BaseElementaryMatrixPtr;

/**
 * @class ElementaryMatrix
 * @brief Class for sd_matr_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryMatrix : public BaseElementaryMatrix {
  private:
    /** @brief Vectors of RESUELEM */
    std::vector< ElementaryTermRealPtr > _realVector;
    std::vector< ElementaryTermComplexPtr > _complexVector;

  public:
    /** @typedef ElementaryMatrixPtr */
    typedef boost::shared_ptr< ElementaryMatrix< ValueType, PhysicalQuantity > >
        ElementaryMatrixPtr;

    /** @brief Constructor with a name */
    ElementaryMatrix( const std::string name )
        : BaseElementaryMatrix(
              name, "MATR_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ){};

    /** @brief Constructor with automatic name */
    ElementaryMatrix() : ElementaryMatrix( ResultNaming::getNewResultName() ){};

    /**
     * @brief Function to update ElementaryTerm
     */
    bool build() {
        if ( _listOfElementaryTerms->exists() ) {
            CALL_JEMARQ();
            _listOfElementaryTerms->updateValuePointer();
            _realVector.clear();
            for ( ASTERINTEGER pos = 0; pos < _listOfElementaryTerms->size(); ++pos ) {
                const std::string name = ( *_listOfElementaryTerms )[pos].toString();
                if ( trim( name ) != "" ) {
                    ElementaryTermRealPtr toPush(
                        boost::make_shared< ElementaryTerm< ASTERDOUBLE > >( name ) );
                    _realVector.push_back( toPush );
                }
            }
            CALL_JEDEMA();
        }
        return true;
    };
};

/** @typedef Elementary matrix for displacement-double */
template class ElementaryMatrix< ASTERDOUBLE, Displacement >;
typedef ElementaryMatrix< ASTERDOUBLE, Displacement > ElementaryMatrixDisplacementReal;
typedef boost::shared_ptr< ElementaryMatrixDisplacementReal > ElementaryMatrixDisplacementRealPtr;

/** @typedef Elementary matrix for displacement-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Displacement >;
typedef ElementaryMatrix< ASTERCOMPLEX, Displacement > ElementaryMatrixDisplacementComplex;
typedef boost::shared_ptr< ElementaryMatrixDisplacementComplex >
    ElementaryMatrixDisplacementComplexPtr;

/** @typedef Elementary matrix for temperature-double */
template class ElementaryMatrix< ASTERDOUBLE, Temperature >;
typedef ElementaryMatrix< ASTERDOUBLE, Temperature > ElementaryMatrixTemperatureReal;
typedef boost::shared_ptr< ElementaryMatrixTemperatureReal > ElementaryMatrixTemperatureRealPtr;

/** @typedef Elementary matrix for pressure-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Pressure >;
typedef ElementaryMatrix< ASTERCOMPLEX, Pressure > ElementaryMatrixPressureComplex;
typedef boost::shared_ptr< ElementaryMatrixPressureComplex > ElementaryMatrixPressureComplexPtr;

#endif /* ELEMENTARYMATRIX_H_ */
