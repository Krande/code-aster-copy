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
#include "MemoryManager/JeveuxUtils.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/DOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryVector
 * @brief Generic class for sd_vect_elem
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
          _listOfLoads( new ListOfLoads() ){};

    /** @brief Constructor with automatic name */
    BaseElementaryVector() : BaseElementaryVector( ResultNaming::getNewResultName() ){};

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

    /** @brief Get the mesh */
    BaseMeshPtr getMesh( void ) const {
        if ( _model ) {
            return _model->getMesh();
        }

        if ( _elemChara ) {
            return _elemChara->getMesh();
        }

        if ( _materialField ) {
            return _materialField->getMesh();
        }

        if ( _listOfLoads ) {
            return _listOfLoads->getMesh();
        }

        return nullptr;
    };

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
typedef std::shared_ptr< BaseElementaryVector > BaseElementaryVectorPtr;

/**
 * @class ElementaryVector
 * @brief Class for sd_vect_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryVector : public BaseElementaryVector {
  private:
    /** @brief Vectors of RESUELEM */
    std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > _elemTerm;

    /** @typedef FieldOnNodesPtr */
    typedef std::shared_ptr< FieldOnNodes< ValueType > > FieldOnNodesPtr;
    FieldOnNodesPtr _veass;

  public:
    /** @typedef ElementaryVectorPtr */
    typedef std::shared_ptr< ElementaryVector< ValueType, PhysicalQuantity > > ElementaryVectorPtr;

    /** @brief Constructor with a name */
    ElementaryVector( const std::string name )
        : BaseElementaryVector(
              name, "VECT_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ){};

    /** @brief Constructor with automatic name */
    ElementaryVector() : ElementaryVector( ResultNaming::getNewResultName() ){};

    /**
     * @brief Function to update ElementaryTerm
     */
    bool build() {
        if ( _elemComp->hasElementaryTerm() ) {

            SetString elemSave;
            for ( auto &elemTerm : _elemTerm ) {
                elemSave.insert( trim( elemTerm->getName() ) );
            }

            auto elemTermNames = _elemComp->getNameOfElementaryTerms();
            SetString elemKeep;
            for ( auto &elemTerm : elemTermNames ) {
                const std::string name = trim( elemTerm.toString() );
                elemKeep.insert( name );
                if ( name != " " && elemSave.count( name ) == 0 ) {
                    std::string name2( name );
                    name2.resize( 19, ' ' );
                    if ( jeveuxExists( name2 + ".REFE" ) ) {
                        // cham_no if .REFE is present, not a resuelem, store in _veass
                        if ( _veass == nullptr or _veass->getName() != name2 ) {
                            _veass = std::make_shared< FieldOnNodes< ValueType > >( name );
                        }
                    } else
                        _elemTerm.push_back(
                            std::make_shared< ElementaryTerm< ValueType > >( name ) );
                }
            }

            // clean ElementaryTerm
            std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > elemTermNew;
            elemTermNew.reserve( _elemTerm.size() );
            for ( auto &elemTerm : _elemTerm ) {
                auto name = trim( elemTerm->getName() );
                if ( elemKeep.count( name ) > 0 ) {
                    elemTermNew.push_back( elemTerm );
                }
            }
            _elemTerm = std::move( elemTermNew );

            if ( !getMesh()->isParallel() && !isMPIFull() ) {
                std::string type = "VECT_ELEM";
                CALLO_SDMPIC( type, getName() );
            }
        }
        _isEmpty = false;
        return true;
    };

    FieldOnNodesPtr getVeass() { return _veass; }

    /**
     * @brief is MPI_COMPLET ?
     */
    bool isMPIFull() {
#ifdef ASTER_HAVE_MPI
        for ( auto &elemTerm : _elemTerm ) {
            if ( !elemTerm->isEmpty() && !elemTerm->isMPIFull() ) {
                return false;
            }
        }
#endif

        return true;
    };

    /**
     * @brief Add elementary term
     */
    void addElementaryTerm( const std::shared_ptr< ElementaryTerm< ValueType > > &elemTerm ) {
        _elemComp->addElementaryTerm( elemTerm->getName() );
        _elemTerm.push_back( elemTerm );
    };

    /**
     * @brief Assembly with dofNume
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesPtr assemble( const BaseDOFNumberingPtr dofNume ) const {
        if ( _isEmpty )
            raiseAsterError( "The ElementaryVector is empty. Call build before" );

        if ( ( !dofNume ) || dofNume->isEmpty() )
            raiseAsterError( "Numerotation is empty" );

        // Create field
        auto field = std::make_shared< FieldOnNodes< ValueType > >( dofNume );

        // Elementary vector names
        std::string vectElemName = getName();
        VectorString vectElemVect( 1, vectElemName );
        char *tabNames = vectorStringAsFStrArray( vectElemVect, 19 );

        // Assembling
        ASTERDOUBLE list_coef = 1.0;
        ASTERINTEGER typscal = typeid( ValueType ) == typeid( ASTERDOUBLE ) ? 1 : 2;
        ASTERINTEGER nbElem = 1;
        std::string base( "G" );

        CALL_ASSVEC( base.c_str(), field->getName().c_str(), &nbElem, tabNames, &list_coef,
                     dofNume->getName().c_str(), &typscal );

        FreeStr( tabNames );

        return field;
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

#endif /* ELEMENTARYVECTOR_H_ */
