#ifndef JEVEUXCOLLECTION_H_
#define JEVEUXCOLLECTION_H_

/**
 * @file JeveuxCollection.h
 * @brief Fichier entete de la classe JeveuxCollection
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "aster_fort_jeveux.h"
#include "aster_utils.h"
#include "astercxx.h"

#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxCollectionObject.h"
#include "MemoryManager/JeveuxObject.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/NamesMap.h"
#include "Utilities/Tools.h"

/**
 * @enum JeveuxCollectionAccessType
 */
enum JeveuxCollectionAccessType { Numbered, Named };

/**
 * @enum JeveuxCollectionMemoryStorageType
 */
enum JeveuxCollectionMemoryStorageType { Sparse, Contiguous };

/**
 * @enum JeveuxCollectionObjectSizes
 */
enum JeveuxCollectionObjectSizes { Constant, Variable };

/**
 * @struct AllowedAccessType
 * @brief Template struct to limit types of access in JeveuxCollectionClass
 * @tparam T Type autorise
 */
template < typename T >
struct AllowedAccessType; // undefined for bad types!

template <>
struct AllowedAccessType< ASTERINTEGER > {
    typedef ASTERINTEGER type;
};

template <>
struct AllowedAccessType< NamesMapChar24 > {
    typedef NamesMapChar24 type;
};

template < typename T, typename U >
using IsSame = typename std::enable_if< std::is_same< T, U >::value >::type;

template < typename T, typename U >
using IsNotSame = typename std::enable_if< !std::is_same< T, U >::value >::type;

/**
 * @class JeveuxCollectionClass
 * @brief Cette classe template permet de definir une collection Jeveux
 * @author Nicolas Sellenet
 */
template < class ValueType, class AccessType = ASTERINTEGER >
class JeveuxCollectionClass : public JeveuxObjectClass, private AllowedAccessType< AccessType > {
  private:
    /** @brief Definition d'un objet de collection du type ValueType */
    typedef JeveuxCollectionObject< ValueType > JeveuxCollObjValType;
    /** @brief std::map associant une chaine a un JeveuxCollObjValType */
    typedef std::map< std::string, ASTERINTEGER > mapStrInt;

    /** @brief La collection est-elle vide ? */
    bool _isEmpty;
    /** @brief La collection est-elle nommée ? */
    bool _isNamed;
    /** @brief Listes de objets de collection */
    std::vector< JeveuxCollObjValType > _listObjects;
    /** @brief Correspondance nom/numéro */
    mapStrInt _mapNumObject;
    ASTERINTEGER _size, _capacity;
    /**
     * @brief Pointeur vers un NamesMap
     * @todo ASTERINTEGER par defaut : pas terrible
     */
    AccessType _namesMap;

    /**
     * @brief Allocation
     */
    bool genericAllocation( ASTERINTEGER size, JeveuxCollectionAccessType access = Named,
                            JeveuxCollectionMemoryStorageType storage = Sparse,
                            JeveuxCollectionObjectSizes objectSizes = Variable,
                            const std::string &name = "", ASTERINTEGER totalSize = 0 ) {
        ASTERINTEGER taille = size;
        _capacity = size;
        _size = 0;
        const std::string strJeveuxBase = "G";
        const ASTERINTEGER intType = AllowedJeveuxType< ValueType >::numTypeJeveux;
        std::string carac( strJeveuxBase + " V " + JeveuxTypesNames[intType] );

        std::string typeColl1( "NO" );
        if ( access == Numbered ) {
            typeColl1 = "NU";
        } else {
            _isNamed = true;
            typeColl1 = typeColl1 + ' ' + name;
        }

        std::string typeColl2( "DISPERSE" );
        if ( storage == Contiguous )
            typeColl2 = "CONTIG";
        std::string typeColl3( "VARIABLE" );
        if ( objectSizes == Constant )
            typeColl3 = "CONSTANT";

        CALLO_JECREC( _name, carac, typeColl1, typeColl2, typeColl3, &taille );
        _isEmpty = false;
        if ( storage == Contiguous ) {
            if ( totalSize <= 0 ) {
                AS_ABORT( "Total size of a contiguous collection must be grower than 0" );
            }
            std::string strParam( "LONT" );
            taille = totalSize;
            CALLO_JEECRA_WRAP( _name, strParam, &taille );
        }

        _listObjects.reserve( _capacity );

        return true;
    };

  public:
    /**
     * @brief Constructeur dans le cas où AccessType n'a pas d'importance
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    JeveuxCollectionClass( const std::string &name )
        : JeveuxObjectClass( name ),
          _isNamed( false ),
          _isEmpty( true ),
          _namesMap( 0 ),
          _size( 0 ),
          _capacity( 0 ) {}

    /**
     * @brief Constructeur dans le cas où AccessType est un NamesMap
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    JeveuxCollectionClass( const std::string &name, AccessType ptr )
        : JeveuxObjectClass( name ),
          _isNamed( false ),
          _isEmpty( true ),
          _namesMap( ptr ),
          _size( 0 ),
          _capacity( 0 ){};

    ~JeveuxCollectionClass(){};

    struct const_iterator {
        ASTERINTEGER position;
        const std::vector< JeveuxCollObjValType > &values;

        inline const_iterator( ASTERINTEGER memoryPosition,
                               const std::vector< JeveuxCollObjValType > &vec )
            : position( memoryPosition ), values( vec ){};

        inline const_iterator( const const_iterator &iter )
            : position( iter.position ), values( iter.values ){};

        inline const_iterator &operator=( const const_iterator &testIter ) {
            position = testIter.position;
            values = testIter.values;
            return *this;
        };

        inline const_iterator &operator++() {
            ++position;
            return *this;
        };

        inline bool operator==( const const_iterator &testIter ) const {
            return testIter.position == position;
        };

        inline bool operator!=( const const_iterator &testIter ) const {
            return !( testIter.position == position );
        };

        inline const JeveuxCollObjValType &operator->() const { return values[position]; };

        inline const JeveuxCollObjValType &operator*() const { return values[position]; };
    };

    /**
     * @brief
     */
    auto begin() const { return _listObjects.begin(); };

    auto end() const { return _listObjects.end(); };

    auto cbegin() const { return _listObjects.cbegin(); };

    auto cend() const { return _listObjects.cend(); };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param access type of access
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    typename std::enable_if< !std::is_same< T1, ASTERINTEGER >::value, bool >::type
    allocateSparseNamed( ASTERINTEGER size, JeveuxCollectionObjectSizes objectSizes = Variable ) {
        if ( !_namesMap->exists() )
            _namesMap->allocate( size );

        if ( _namesMap->capacity() != size ) {
            AS_ABORT( "Sizes do not match: " + std::to_string( size ) + " vs " +
                      std::to_string( _namesMap->capacity() ) );
        }

        return genericAllocation( size, Named, Sparse, objectSizes, _namesMap->getName() );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param access type of access
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value, bool >::type
    allocateSparseNumbered( ASTERINTEGER size,
                            JeveuxCollectionObjectSizes objectSizes = Variable ) {
        return genericAllocation( size, Numbered, Sparse, objectSizes, "" );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param totalSize total size of the collection
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    typename std::enable_if< !std::is_same< T1, ASTERINTEGER >::value, bool >::type
    allocateContiguousNamed( ASTERINTEGER size, ASTERINTEGER totalSize,
                             JeveuxCollectionObjectSizes objectSizes = Variable ) {
        if ( !_namesMap->exists() )
            _namesMap->allocate( size );

        if ( _namesMap->capacity() != size ) {
            AS_ABORT( "Sizes do not match: " + std::to_string( size ) + " vs " +
                      std::to_string( _namesMap->capacity() ) );
        }

        return genericAllocation( size, Named, Contiguous, objectSizes, _namesMap->getName(),
                                  totalSize );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param totalSize total size of the collection
     * @param access type of access
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value, bool >::type
    allocateContiguousNumbered( ASTERINTEGER size, ASTERINTEGER totalSize,
                                JeveuxCollectionObjectSizes objectSizes = Variable ) {
        return genericAllocation( size, Numbered, Contiguous, objectSizes, "", totalSize );
    };

    /**
     * @brief Allocation of one object by name
     */
    JeveuxCollObjValType allocateObject( const std::string &name, const ASTERINTEGER &nbValues ) {
        if ( size() >= capacity() ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( size() ) + " vs " +
                      std::to_string( capacity() ) );
        }

        if ( !isNamed() ) {
            AS_ABORT( "Collection is not named" );
        }

        _mapNumObject[std::string( trim( name.c_str() ) )] = size();

        _size++;
        JeveuxCollObjValType obj( _name, size(), name, nbValues );

        _listObjects.push_back( obj );
        return obj;
    };

    /**
     * @brief Allocation of one object by name
     */
    JeveuxCollObjValType allocateObject( const std::string &name,
                                         const std::vector< ValueType > &values ) {
        auto obj = allocateObject( name, values.size() );
        obj.setValues( values );

        return obj;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    JeveuxCollObjValType allocateObject( const ASTERINTEGER &nbValues ) {
        if ( size() >= capacity() ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( size() ) + " vs " +
                      std::to_string( capacity() ) );
        }

        if ( !isNumeroted() ) {
            AS_ABORT( "Collection is not numeroted" );
        }

        _size++;
        JeveuxCollObjValType obj( _name, size(), nbValues );

        _listObjects.push_back( obj );
        return obj;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    JeveuxCollObjValType allocateObject( const std::vector< ValueType > &values ) {
        auto obj = allocateObject( values.size() );
        obj.setValues( values );

        return obj;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    void push_back( const std::vector< ValueType > &values ) {
        auto obj = allocateObject( values.size() );
        obj.setValues( values );
    };
    /**
     * @brief Allocation of one object by name
     */
    void push_back( const std::string &name, const std::vector< ValueType > &values ) {
        auto obj = allocateObject( name, values.size() );
        obj.setValues( values );
    };

    /**
     * @brief Deallocate collection
     */
    void deallocate() {
        if ( _name != "" && get_sh_jeveux_status() == 1 )
            CALLO_JEDETR( _name );

        _size = 0;
        _capacity = 0;
        _listObjects.clear();
    };

    /**
     * @brief Methode permettant de construire une collection a partir d'une collection
     *   existante en memoire Jeveux
     * @return Renvoit true si la construction s'est bien deroulee
     */
    bool build( bool force = false );

    /**
     * @brief Methode verifiant l'existence d'un objet de collection dans la collection
     * @param name Chaine contenant le nom de l'objet
     * @return Renvoit true si l'objet existe dans la collection
     */
    bool existsObject( const std::string &name ) const;

    /**
     * @brief Methode verifiant l'existence d'un objet de collection dans la collection
     * @param number entier
     * @return Renvoit true si l'objet existe dans la collection
     */
    bool existsObject( const ASTERINTEGER &number ) const;

    /**
     * @brief Methode permettant d'obtenir la liste des objets nommés dans la collection
     * @return vecteur de noms d'objet de collection
     */
    std::vector< JeveuxChar32 > getObjectsNames() const;

    const std::vector< JeveuxCollObjValType > &getObjects() const { return _listObjects; };

    inline const JeveuxCollObjValType &operator[]( const ASTERINTEGER &position ) const {
#ifdef ASTER_DEBUG_CXX
        if ( _isEmpty ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( position > size() || position <= 0 ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( position ) + " vs " +
                      std::to_string( size() ) );
        }
#endif
        return _listObjects[position - 1];
    };

    inline JeveuxCollObjValType &operator[]( const ASTERINTEGER &position ) {
#ifdef ASTER_DEBUG_CXX
        if ( _isEmpty ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( position > size() || position <= 0 ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( position ) + " vs " +
                      std::to_string( size() ) );
        }
#endif
        return _listObjects[position - 1];
    };

    inline const JeveuxCollObjValType &operator[]( const std::string &name ) const {
#ifdef ASTER_DEBUG_CXX
        if ( _isEmpty ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( !_isNamed ) {
            AS_ABORT( "Collection " + _name + " is not named" );
        }
#endif
        const auto &curIter = _mapNumObject.find( trim( name ) );
        if ( curIter == _mapNumObject.end() ) {
            AS_ABORT( "Name not in collection: " + name );
        }

        return _listObjects[curIter->second];
    };

    inline JeveuxCollObjValType &operator[]( const std::string &name ) {
#ifdef ASTER_DEBUG_CXX
        if ( _isEmpty ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( !_isNamed ) {
            AS_ABORT( "Collection " + _name + " is not named" );
        }
#endif
        const auto &curIter = _mapNumObject.find( trim( name ) );
        if ( curIter == _mapNumObject.end() ) {
            AS_ABORT( "Name not in collection: " + name );
        }

        return _listObjects[curIter->second];
    };

    inline ASTERINTEGER size() const { return _size; };

    inline ASTERINTEGER capacity() const { return _capacity; };

    inline ASTERBOOL isNamed() const { return _isNamed; };

    inline ASTERBOOL isNumeroted() const { return ( !_isNamed ); };

    /**
     * @brief Surcharge de l'operateur =
     */
    JeveuxCollectionClass &operator=( JeveuxCollectionClass< ValueType, AccessType > &toCopy ) {

        std::string base( "G" );
        bool dupcol = true;
        CALLO_JEDUPO( toCopy.getName(), base, getName(), (ASTERLOGICAL *)&dupcol );
        AS_ASSERT( build() );

        return *this;
    };

    bool operator==( JeveuxCollectionClass< ValueType, AccessType > &toCompar ) {

        if ( this->size() != toCompar.size() )
            return false;

        auto size = this->size();

        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            if ( _listObjects[i] == toCompar._listObjects[i] ) {
                // nothing
            } else {
                return false;
            }
        }

        return true;
    };

    JeveuxCollectionClass< ValueType, AccessType > &operator*=( const ValueType &scal ) {
        auto size = this->size();

        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            _listObjects[i] *= scal;
        }

        return *this;
    };
};

template < class ValueType, class AccessType >
bool JeveuxCollectionClass< ValueType, AccessType >::build( bool force ) {

    if ( !exists() ) {
        return false;
    }

    ASTERINTEGER nbColObj, valTmp;
    JeveuxChar8 param( "NMAXOC" );
    std::string charval( 32, ' ' );
    CALLO_JELIRA( _name, param, &nbColObj, charval );

    if ( !force && !_isEmpty && size() == nbColObj ) {
        for ( ASTERINTEGER i = 0; i < nbColObj; ++i )
            if ( _listObjects[i].hasMoved() ) {
                force = true;
                break;
            }
        if ( !force )
            return true;
    }

    _size = nbColObj;
    _capacity = nbColObj;
    _listObjects.clear();
    _listObjects.reserve( _capacity );

    param = "ACCES ";
    charval = std::string( 32, ' ' );
    CALLO_JELIRA( _name, param, &valTmp, charval );
    const std::string resu( charval, 0, 2 );

    if ( resu == "NO" ) {
        _isNamed = true;
    }

    for ( ASTERINTEGER i = 1; i <= _size; ++i ) {
        if ( !existsObject( i ) ) {
            _listObjects.push_back( JeveuxCollObjValType( _name, i, _isNamed, false ) );
        } else {
            _listObjects.push_back( JeveuxCollObjValType( _name, i, _isNamed, true ) );
        }

        if ( _isNamed ) {
            _mapNumObject[_listObjects[i - 1].getStringName()] = i - 1;
        }
    }
    _isEmpty = false;
    return true;
};

template < class ValueType, class AccessType >
bool JeveuxCollectionClass< ValueType, AccessType >::existsObject( const std::string &name ) const {
    std::string charJeveuxName( 32, ' ' );
    ASTERINTEGER returnBool;
    CALLO_JEXNOM( charJeveuxName, _name, name );
    CALLO_JEEXIN( charJeveuxName, &returnBool );
    if ( returnBool == 0 )
        return false;
    return true;
};

template < class ValueType, class AccessType >
bool JeveuxCollectionClass< ValueType, AccessType >::existsObject(
    const ASTERINTEGER &number ) const {
    std::string charJeveuxName( 32, ' ' );
    ASTERINTEGER returnBool = number;
    CALLO_JEXNUM( charJeveuxName, _name, &returnBool );
    CALLO_JEEXIN( charJeveuxName, &returnBool );

    if ( returnBool == 0 )
        return false;
    return true;
};

template < class ValueType, class AccessType >
std::vector< JeveuxChar32 >
JeveuxCollectionClass< ValueType, AccessType >::getObjectsNames() const {
    std::vector< JeveuxChar32 > toReturn;
    toReturn.reserve( size() );
    for ( auto curObject : _listObjects )
        toReturn.push_back( curObject.getName() );
    return toReturn;
};

/**
 * @class JeveuxCollection
 * @brief Enveloppe d'un pointeur intelligent vers un JeveuxCollectionClass
 * @author Nicolas Sellenet
 */
template < class ValueType, class AccessType = ASTERINTEGER >
class JeveuxCollection {
  public:
    typedef std::shared_ptr< JeveuxCollectionClass< ValueType, AccessType > >
        JeveuxCollectionTypePtr;

  private:
    JeveuxCollectionTypePtr _jeveuxCollectionPtr;

  public:
    JeveuxCollection() : _jeveuxCollectionPtr( nullptr ) {}

    /**
     * @brief Constructeur dans le cas où AccessType n'a pas d'importance
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    JeveuxCollection( const std::string &nom )
        : _jeveuxCollectionPtr(
              std::make_shared< JeveuxCollectionClass< ValueType, AccessType > >( nom ) ) {}

    /**
     * @brief Constructeur dans le cas où AccessType est un NamesMap
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    JeveuxCollection( const std::string &nom, AccessType ptr )
        : _jeveuxCollectionPtr(
              std::make_shared< JeveuxCollectionClass< ValueType, AccessType > >( nom, ptr ) ){};

    ~JeveuxCollection(){};

    JeveuxCollection &operator=( const JeveuxCollection< ValueType, AccessType > &tmp ) {
        _jeveuxCollectionPtr = tmp._jeveuxCollectionPtr;
        return *this;
    };

    const JeveuxCollectionTypePtr &operator->() const { return _jeveuxCollectionPtr; };

    JeveuxCollectionClass< ValueType, AccessType > &operator*( void ) const {
        return *_jeveuxCollectionPtr;
    };

    bool isEmpty() const {
        if ( _jeveuxCollectionPtr.use_count() == 0 )
            return true;
        return false;
    };

    auto begin() { return _jeveuxCollectionPtr->begin(); };

    auto end() { return _jeveuxCollectionPtr->end(); };

    auto cbegin() { return _jeveuxCollectionPtr->cbegin(); };

    auto cend() { return _jeveuxCollectionPtr->cend(); };
};

/** @typedef Definition d'une collection de type entier long */
typedef JeveuxCollection< ASTERINTEGER > JeveuxCollectionLong;
/** @typedef Definition d'une collection de type entier court */
typedef JeveuxCollection< ASTERINTEGER4 > JeveuxCollectionShort;
/** @typedef Definition d'une collection de type double */
typedef JeveuxCollection< ASTERDOUBLE > JeveuxCollectionReal;
/** @typedef Definition d'une collection de type double complex */
typedef JeveuxCollection< ASTERCOMPLEX > JeveuxCollectionComplex;
/** @typedef Definition d'une collection de JeveuxChar8 */
typedef JeveuxCollection< JeveuxChar8 > JeveuxCollectionChar8;
/** @typedef Definition d'une collection de JeveuxChar16 */
typedef JeveuxCollection< JeveuxChar16 > JeveuxCollectionChar16;
/** @typedef Definition d'une collection de JeveuxChar24 */
typedef JeveuxCollection< JeveuxChar24 > JeveuxCollectionChar24;
/** @typedef Definition d'une collection de JeveuxChar32 */
typedef JeveuxCollection< JeveuxChar32 > JeveuxCollectionChar32;
/** @typedef Definition d'une collection de JeveuxChar80 */
typedef JeveuxCollection< JeveuxChar80 > JeveuxCollectionChar80;

#endif /* JEVEUXCOLLECTION_H_ */
