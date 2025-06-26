#ifndef NAMESMAP_H_
#define NAMESMAP_H_

/**
 * @file NamesMap.h
 * @brief Fichier entete de la classe JeveuxBidirectionnalMap
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "aster_fort_jeveux.h"
#include "aster_utils.h"

#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxObject.h"

/**
 * @class NamesMapClass
 * @brief Equivalent du pointeur de nom dans Jeveux
 * @author Nicolas Sellenet
 */
template < typename ValueType >
class NamesMapClass : public JeveuxObjectClass, private AllowedJeveuxType< ValueType > {
  private:
  public:
    /**
     * @brief Constructeur
     * @param name Nom Jeveux de l'objet
     */
    NamesMapClass( std::string name ) : JeveuxObjectClass( name ) {};

    /**
     * @brief Destructeur
     */
    ~NamesMapClass() {};

    /**
     * @brief Ajout d'un élément
     * @param position position of element to add
     * @param toAdd value to add
     * @return true if adding is ok
     */
    bool add( const ASTERINTEGER &position, const ValueType &toAdd ) {
        if ( position <= this->capacity() ) {
            JeveuxChar32 objName( " " );
            CALLO_JEXNOM( objName, _name, toAdd );
            CALLO_JECROC( objName );
            return true;
        } else {
            std::string error = "Out of range of NamesMap '" + _name +
                                "', index = " + std::to_string( position ) +
                                " ( capacity = " + std::to_string( this->capacity() ) + " )";
            throw std::out_of_range( error );
        }
        return false;
    };

    /**
     * @brief Allocation
     * @param size Taille
     * @return vrai en cas d'allocation
     */
    void allocate( ASTERINTEGER size ) {
        if ( _name != "" && size > 0 ) {
            std::string strJeveuxBase = JeveuxMemoryTypesNames[_mem];
            ASTERINTEGER taille = size;
            const auto intType = AllowedJeveuxType< ValueType >::numTypeJeveux;
            std::string carac = strJeveuxBase + " N " + JeveuxTypesNames[intType];
            CALLO_JECREO( _name, carac );
            std::string param( "NOMMAX" );
            CALLO_JEECRA_WRAP( _name, param, &taille );
        } else {
            throw std::bad_alloc();
        }
    };

    /**
     * @brief Desallocation d'un vecteur Jeveux
     */
    void deallocate() {
        if ( _name != "" && get_sh_jeveux_status() == 1 )
            CALLO_JEDETR( _name );
    };

    /**
     * @brief Recuperation de la chaine correspondante a l'entier
     * @param index Numero de l'element demande
     * @return Chaine de caractere correspondante
     */
    std::string getStringFromIndex( ASTERINTEGER index ) const {
        JeveuxChar32 objName( " " );
        JeveuxChar32 charName( " " );
        CALLO_JEXNUM( objName, _name, &index );
        CALLO_JENUNO( objName, charName );
        return charName.toString();
    };

    /**
     * @brief Recuperation de l'entier correspondant a une chaine
     * @param string Chaine recherchee
     * @return Entier correspondant
     */
    ASTERINTEGER getIndexFromString( const std::string &string ) const {
        JeveuxChar32 objName( " " );
        CALLO_JEXNOM( objName, _name, string );
        ASTERINTEGER resu = -1;
        CALLO_JENONU( objName, &resu );
        return resu;
    };

    /**
     * @brief Get the size (=NOMUTI)
     * @return size of object
     */
    ASTERINTEGER size() const {
        if ( !exists() )
            return 0;

        ASTERINTEGER vectSize;
        JeveuxChar8 param( "NOMUTI" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( _name, param, &vectSize, dummy );
        return vectSize;
    };

    /**
     * @brief Get the capacity  (=NOMMAX)
     * @return capicity of object
     */
    ASTERINTEGER capacity() const {
        if ( !exists() )
            return 0;

        ASTERINTEGER vectSize;
        JeveuxChar8 param( "NOMMAX" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( _name, param, &vectSize, dummy );
        return vectSize;
    };

    bool operator==( NamesMapClass< ValueType > &toCompare ) {
        if ( this->size() != toCompare.size() )
            return false;

        const auto size = toCompare.size();
        for ( ASTERINTEGER i = 1; i <= size; ++i ) {
            if ( this->getStringFromIndex( i ) != toCompare.getStringFromIndex( i ) )
                return false;
        }
        return true;
    };
};

/**
 * class NamesMap
 *   Enveloppe d'un pointeur intelligent vers un NamesMapClass
 * @author Nicolas Sellenet
 */
template < class ValueType >
class NamesMap {
  public:
    typedef std::shared_ptr< NamesMapClass< ValueType > > NamesMapPtr;

  private:
    NamesMapPtr _namesMapPtr;

  public:
    NamesMap( std::string nom ) : _namesMapPtr( new NamesMapClass< ValueType >( nom ) ) {};

    ~NamesMap() {};

    NamesMap &operator=( const NamesMap< ValueType > &tmp ) {
        _namesMapPtr = tmp._namesMapPtr;
        return *this;
    };

    const NamesMapPtr &operator->( void ) const { return _namesMapPtr; };

    NamesMapClass< ValueType > &operator*( void ) const { return *_namesMapPtr; };

    bool exists() const {
        if ( _namesMapPtr == nullptr )
            return false;

        if ( _namesMapPtr.use_count() == 0 )
            return false;

        return _namesMapPtr->exists();
    }
};

/** @typedef Definition d'un pointeur de nom Jeveux entier long */
typedef NamesMap< ASTERINTEGER > NamesMapLong;
/** @typedef Definition d'un pointeur de nom Jeveux entier court */
typedef NamesMap< ASTERINTEGER4 > NamesMapShort;
/** @typedef Definition d'un pointeur de nom Jeveux double */
typedef NamesMap< ASTERDOUBLE > NamesMapReal;
/** @typedef Definition d'un pointeur de nom Jeveux double complex */
typedef NamesMap< ASTERCOMPLEX > NamesMapComplex;
/** @typedef Definition d'un vecteur de JeveuxChar8 */
typedef NamesMap< JeveuxChar8 > NamesMapChar8;
/** @typedef Definition d'un pointeur de nom JeveuxChar16 */
typedef NamesMap< JeveuxChar16 > NamesMapChar16;
/** @typedef Definition d'un pointeur de nom JeveuxChar24 */
typedef NamesMap< JeveuxChar24 > NamesMapChar24;
/** @typedef Definition d'un pointeur de nom JeveuxChar32 */
typedef NamesMap< JeveuxChar32 > NamesMapChar32;
/** @typedef Definition d'un pointeur de nom JeveuxChar80 */
typedef NamesMap< JeveuxChar80 > NamesMapChar80;
/** @typedef Definition d'un pointeur de nom JeveuxLogical */
typedef NamesMap< ASTERBOOL > NamesMapLogical;

#endif /* NAMESMAP_H_ */
