#ifndef ELEMENTARYTERM_H_
#define ELEMENTARYTERM_H_

/**
 * @file ElementaryTerm.h
 * @brief Fichier entete de la classe ElementaryTerm
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "aster_fort_utils.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"

/**
 * @class ElementaryTerm
 * @brief Class which describe a RESUELEM (which is part of MATR_ELEM and VECT_ELEM )
 */
template < typename ValueType >
class ElementaryTerm : public DataField {
  private:
    /** @brief Objet Jeveux '.NOLI' */
    JeveuxVectorChar24 _noli;
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Objet Jeveux '.RESL' */
    JeveuxCollection< ValueType > _resl;
    /** @brief Ligel */
    FiniteElementDescriptorPtr _FEDesc;

  public:
    /**
     * @brief Constructor only with predefined name
     * @param name predefined name
     */
    ElementaryTerm( const std::string name )
        : DataField( name, "RESUELEM" ),
          _noli( JeveuxVectorChar24( getName() + ".NOLI" ) ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _resl( JeveuxCollection< ValueType >( getName() + ".RESL" ) ){};

    ElementaryTerm() : ElementaryTerm( DataStructureNaming::getNewName() ){};

    void setFiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc ) {
        _FEDesc = FEDesc;
    };

    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    bool exists() const { return _noli->exists() && _descriptor->exists() && _resl->exists(); };

    std::string getOption() {
        _noli->updateValuePointer();
        return trim( ( *_noli )[1].toString() );
    };

    BaseMeshPtr getMesh() const {
        if ( _FEDesc ) {
            return _FEDesc->getMesh();
        }
        return nullptr;
    }

    std::string getPhysicalQuantity() {
        const std::string typeco( "RESUELEM" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "NOM_GD" );
        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
        return trim( repk.toString() );
    }

    ASTERINTEGER getNumberOfGroupOfCells() const {
        _descriptor->updateValuePointer();

        return ( *_descriptor )[1];
    }

    /**
     * @brief Return MODE_LOCAL
     */
    std::string getLocalMode() const {
        const auto nbGrel = getNumberOfGroupOfCells();
        _descriptor->updateValuePointer();
        const std::string cata = "&CATA.TE.NOMMOLOC";
        JeveuxChar24 objName, charName;

        std::string modeName;
        for ( auto igr = 0; igr < nbGrel; igr++ ) {
            auto mode = ( *_descriptor )[2 + igr];
            if ( mode > 0 ) {
                CALLO_JEXNUM( objName, cata, &mode );
                CALLO_JENUNO( objName, charName );
                auto modeN = trim( charName.toString().substr( 14, 10 ) );
                if ( modeName.empty() ) {
                    modeName = modeN;
                } else {
                    if ( modeName != modeN ) {
                        AS_ABORT( "Multiple names." );
                    }
                }
            }
        }

        return modeName;
    }

    bool isMPIFull() {
        AS_ASSERT( _noli->exists() );
        _noli->updateValuePointer();
        return trim( ( *_noli )[2].toString() ) == "MPI_COMPLET";
    }

    bool isEmpty() { return !( _noli->exists() && _descriptor->exists() ); }

    bool build() { return _resl->build( true ); };

    /**
     * @brief TimesEqual overloading
     */
    ElementaryTerm< ValueType > &operator*=( const ValueType &scal ) {

        ( *_resl ) *= scal;

        return *this;
    };
};

/** @typedef ElementaryTermRealPtr */
using ElementaryTermReal = ElementaryTerm< ASTERDOUBLE >;
using ElementaryTermRealPtr = std::shared_ptr< ElementaryTermReal >;

/** @typedef ElementaryTermComplexPtr */
using ElementaryTermComplex = ElementaryTerm< ASTERCOMPLEX >;
using ElementaryTermComplexPtr = std::shared_ptr< ElementaryTermComplex >;

#endif /* ELEMENTARYTERM_H_ */
