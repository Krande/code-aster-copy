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

#include "DataFields/ElementaryTerm.h"
#include "LinearAlgebra/BaseElementaryVector.h"
#include "MemoryManager/JeveuxUtils.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class GenericElementaryVector
 * @brief Generic class for sd_vect_elem
 */

template < typename ValueType >
class GenericElementaryVector : public BaseElementaryVector {
  protected:
    /** @brief Vectors of RESUELEM */
    std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > _elemTerm;

    /** @typedef FieldOnNodesPtr */
    typedef std::shared_ptr< FieldOnNodes< ValueType > > FieldOnNodesPtr;

    // FIXME: should not exist
    FieldOnNodesPtr _veass;

  public:
    /** @typedef GenericElementaryVectorPtr */
    typedef std::shared_ptr< GenericElementaryVector< ValueType > > GenericElementaryVectorPtr;

    /** @brief Constructor with a name */
    GenericElementaryVector( const std::string name,
                             const std::string type = "VECT_ELEM_" +
                                                      std::string( typeid( ValueType ) ==
                                                                           typeid( ASTERDOUBLE )
                                                                       ? "_R"
                                                                       : "_C" ) )
        : BaseElementaryVector( name, type ){};

    /** @brief Constructor with automatic name */
    GenericElementaryVector() : GenericElementaryVector( ResultNaming::getNewResultName() ){};

    GenericElementaryVector( const PhysicalProblemPtr phys_pb ) : GenericElementaryVector() {
        this->setPhysicalProblem( phys_pb );
    };

    /**
     * @brief Return MODE_LOCAL
     */
    std::string getLocalMode() const {
        if ( _elemTerm.size() == 0 ) {
            return std::string();
        }

        std::string modeName = _elemTerm[0]->getLocalMode();

        for ( auto &term : _elemTerm ) {
            auto modeN = term->getLocalMode();
            if ( modeName != modeN ) {
                AS_ABORT( "Multiple names." );
            }
        }

        return modeName;
    }

    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptor() const {
        std::vector< FiniteElementDescriptorPtr > FEDs;

        FEDs.reserve( _elemTerm.size() );
        SetString elemKeep;

        for ( auto &elemTerm : _elemTerm ) {
            const auto FED = elemTerm->getFiniteElementDescriptor();
            if ( FED ) {
                const auto name = FED->getName();
                if ( name != " " && elemKeep.count( name ) == 0 ) {
                    FEDs.push_back( FED );
                    elemKeep.insert( name );
                }
            }
        }

        // if ( _model ) {
        //     if ( elemKeep.count( _model->getName() ) == 0 ) {
        //         FEDs.push_back( _model->getFiniteElementDescriptor() );
        //         elemKeep.insert( _model->getName() );
        //     }
        // }

        return FEDs;
    }

    /** @brief Get the mesh */
    BaseMeshPtr getMesh( void ) const {

        auto FEDs = getFiniteElementDescriptor();
        for ( auto &FED : FEDs ) {
            if ( FED && FED->getMesh() ) {
                return FED->getMesh();
            }
        }

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
     * @brief Add elementary term
     */
    void addElementaryTerm( const std::shared_ptr< ElementaryTerm< ValueType > > &elemTerm,
                            const ASTERINTEGER iload ) {
        CALLO_CORICHWRITE( elemTerm->getName(), &iload );
        _elemComp->addElementaryTerm( elemTerm->getName() );
        _elemTerm.push_back( elemTerm );
    };

    /**
     * @brief Assembly with dofNume
     * @param dofNume object DOFNumbering
     */
    FieldOnNodesPtr assemble( const BaseDOFNumberingPtr dofNume = nullptr ) const {
        if ( _isEmpty )
            raiseAsterError( "The ElementaryVector is empty. Call build before" );

        BaseDOFNumberingPtr nume = dofNume;

        if ( ( !nume ) || nume->isEmpty() ) {
            nume = std::make_shared< DOFNumbering >();
            nume->addFiniteElementDescriptors( getFiniteElementDescriptor() );
            nume->computeNumberingWithLocalMode( getLocalMode() );
        }

        // Create field
        auto field = std::make_shared< FieldOnNodes< ValueType > >( nume );

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
                     nume->getName().c_str(), &typscal );

        FreeStr( tabNames );

        return field;
    };
};

/** @typedef GenericElementary vector for double */
template class GenericElementaryVector< ASTERDOUBLE >;
using ElementaryVectorReal = GenericElementaryVector< ASTERDOUBLE >;
using ElementaryVectorRealPtr = std::shared_ptr< ElementaryVectorReal >;

/** @typedef GenericElementary vector for complex */
template class GenericElementaryVector< ASTERCOMPLEX >;
using ElementaryVectorComplex = GenericElementaryVector< ASTERCOMPLEX >;
using ElementaryVectorComplexPtr = std::shared_ptr< ElementaryVectorComplex >;
