#ifndef JEVEUXCOLLECTION_H_
#define JEVEUXCOLLECTION_H_

/**
 * @file JeveuxCollection.h
 * @brief Fichier entete de la classe JeveuxCollection
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
#include "MemoryManager/JeveuxCollectionObject.h"
#include "MemoryManager/JeveuxObject.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/NamesMap.h"
#include "Utilities/MapVectorEnumerator.h"
#include "Utilities/Tools.h"

#include <type_traits>

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

/**
 * @struct CollectionMemoryType
 * @brief Template struct to specify collection memory type
 * @tparam T Type autorise
 */
template < JeveuxCollectionMemoryStorageType T >
struct CollectionMemoryType;

template <>
struct CollectionMemoryType< Sparse > {
    typedef ASTERINTEGER type;
};

template <>
struct CollectionMemoryType< Contiguous > {
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
template < class ValueType, class AccessType = ASTERINTEGER,
           JeveuxCollectionMemoryStorageType MemoryType = Sparse >
class JeveuxCollectionClass : public JeveuxObjectClass,
                              private AllowedAccessType< AccessType >,
                              private CollectionMemoryType< MemoryType > {
  private:
    /** @brief Definition d'un objet de collection du type ValueType */
    typedef JeveuxCollectionObject< ValueType > JeveuxCollObjValType;
    /** @brief std::map associant une chaine a un JeveuxCollObjValType */
    typedef std::map< std::string, ASTERINTEGER > mapStrInt;

    /** @brief La collection a-t-elle été construite ? */
    bool _isBuilt;
    /** @brief La collection est-elle nommée ? */
    bool _isNamed;
    /** @brief Nombre d'objets de collection */
    ASTERINTEGER _objectCount = 0;
    /** @brief Correspondance nom/numéro */
    mapStrInt _mapNumObject;
    ASTERINTEGER _capacity;
    ASTERINTEGER _totalSize;
    JeveuxCollObjValType _temporaryCollectionObject;
    JeveuxCollObjValType _temporaryEndCollectionObject;
    VectorLong _indirection;
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
        _isBuilt = true;
        if ( storage == Contiguous ) {
            if ( totalSize <= 0 ) {
                AS_ABORT( "Total size of a contiguous collection must be grower than 0" );
            }
            std::string strParam( "LONT" );
            taille = totalSize;
            CALLO_JEECRA_WRAP( _name, strParam, &taille );
            _totalSize = totalSize;
        }

        return true;
    };

    ASTERINTEGER getMaxIndex() const { return _objectCount; };

    ASTERINTEGER getNewIndex() const {
        auto newIndex = getMaxIndex() + 1;
        if ( newIndex > _capacity ) {
            AS_ABORT( "Can not add new object" );
        }

        return newIndex;
    };

  public:
    // template< class ValueType, class AccessType, JeveuxCollectionMemoryStorageType MemoryType >
    class FastCollectionObject {
        const int _position;
        int _size;
        const JeveuxCollectionClass< ValueType, AccessType, MemoryType > &_ref;
        ValueType *_ptr = nullptr;

      public:
        FastCollectionObject( const JeveuxCollectionClass< ValueType, AccessType, MemoryType > &ref,
                              const int &pos )
            : _ref( ref ), _position( pos ) {
            std::string charJeveuxName( 32, ' ' );
            ASTERINTEGER num = _position;
            CALLO_JEXNUM( charJeveuxName, _ref.getName(), &num );

            const std::string read( "L" );
            CALLO_JEVEUOC( charJeveuxName, read, (void *)( &_ptr ) );

            ASTERINTEGER valTmp;
            JeveuxChar8 param( "LONMAX" );
            std::string charval = std::string( 32, ' ' );
            CALLO_JELIRA( charJeveuxName, param, &valTmp, charval );
            _size = valTmp;
        };

        std::vector< ValueType > toVector() const {
            std::vector< ValueType > toReturn( _size, 0 );
            std::copy( _ptr, _ptr + _size, &toReturn[0] );
            return toReturn;
        };
    };

    /**
     * @brief Constructeur dans le cas où AccessType n'a pas d'importance
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    JeveuxCollectionClass( const std::string &name )
        : JeveuxObjectClass( name ),
          _isNamed( false ),
          _isBuilt( false ),
          _namesMap( 0 ),
          _capacity( 0 ),
          _totalSize( 0 ),
          _temporaryCollectionObject(),
          _temporaryEndCollectionObject() {}

    /**
     * @brief Constructeur dans le cas où AccessType est un NamesMap
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    JeveuxCollectionClass( const std::string &name, AccessType ptr )
        : JeveuxObjectClass( name ),
          _isNamed( false ),
          _isBuilt( false ),
          _namesMap( ptr ),
          _capacity( 0 ),
          _totalSize( 0 ),
          _temporaryCollectionObject(),
          _temporaryEndCollectionObject() {};

    ~JeveuxCollectionClass() {
        // #ifdef ASTER_DEBUG_CXX_OBJECTS
        //         std::cout << "DEBUG: JeveuxCollection.destr: " << _name <<
        //         std::endl;
        // #endif
        this->deallocate();
    };

    /** @brief Iterator on collection */
    template < class ValueTypeIter = ValueType, class AccessTypeIter = AccessType,
               JeveuxCollectionMemoryStorageType MemoryTypeIter = MemoryType >
    class iterator {
        const JeveuxCollectionClass< ValueTypeIter, AccessTypeIter, MemoryTypeIter > *_collection;
        typedef std::pair< int, JeveuxCollectionObject< ValueTypeIter > > iterPair;
        int _index = 1;

      public:
        iterator(
            const JeveuxCollectionClass< ValueTypeIter, AccessTypeIter, MemoryTypeIter > *coll,
            const int &position )
            : _collection( coll ), _index( position ) {}

        bool operator!=( const iterator &other ) const {
            if ( _collection->getName() != other._collection->getName() ) {
                return false;
            } else {
                return _index != other._index;
            }
        }
        void operator++() { ++_index; }
        auto operator*() {
            return iterPair( _index,
                             JeveuxCollectionObject< ValueTypeIter >( ( *_collection )[_index] ) );
        }
    };

    /** @brief Constant iterator on collection */
    template < class ValueTypeIter = ValueType, class AccessTypeIter = AccessType,
               JeveuxCollectionMemoryStorageType MemoryTypeIter = MemoryType >
    class const_iterator {
        const JeveuxCollectionClass< ValueTypeIter, AccessTypeIter, MemoryTypeIter > *_collection;
        typedef std::pair< int, const JeveuxCollectionObject< ValueTypeIter > > iterPair;
        int _index = 1;

      public:
        const_iterator(
            const JeveuxCollectionClass< ValueTypeIter, AccessTypeIter, MemoryTypeIter > *coll,
            const int &position )
            : _collection( coll ), _index( position ) {}

        bool operator!=( const const_iterator &other ) const {
            if ( _collection->getName() != other._collection->getName() ) {
                return false;
            } else {
                return _index != other._index;
            }
        }
        void operator++() { ++_index; }
        const auto &operator*() const {
            return iterPair( _index,
                             JeveuxCollectionObject< ValueTypeIter >( ( *_collection )[_index] ) );
        }
    };

    /**
     * @brief
     */
    auto begin() const {
        const auto firstPos = _indirection[0];
        return iterator< ValueType, AccessType, MemoryType >( this, firstPos );
    };

    auto end() const {
        const auto lastPos = _indirection[_objectCount - 1];
        return iterator< ValueType, AccessType, MemoryType >( this, lastPos + 1 );
    };

    auto cbegin() const {
        const auto firstPos = _indirection[0];
        return const_iterator< ValueType, AccessType, MemoryType >( this, firstPos );
    };

    auto cend() const {
        const auto lastPos = _indirection[_objectCount - 1];
        return const_iterator< ValueType, AccessType, MemoryType >( this, lastPos + 1 );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param access type of access
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    typename std::enable_if< !std::is_same< T1, ASTERINTEGER >::value && MemoryType == Sparse,
                             bool >::type
    allocate( ASTERINTEGER size, JeveuxCollectionObjectSizes objectSizes = Variable ) {
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
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value && MemoryType == Sparse,
                             bool >::type
    allocate( ASTERINTEGER size, JeveuxCollectionObjectSizes objectSizes = Variable ) {
        return genericAllocation( size, Numbered, Sparse, objectSizes, "" );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param totalSize total size of the collection
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    typename std::enable_if< !std::is_same< T1, ASTERINTEGER >::value && MemoryType == Contiguous,
                             bool >::type
    allocate( ASTERINTEGER size, ASTERINTEGER totalSize,
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
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value && MemoryType == Contiguous,
                             bool >::type
    allocate( ASTERINTEGER size, ASTERINTEGER totalSize,
              JeveuxCollectionObjectSizes objectSizes = Variable ) {
        return genericAllocation( size, Numbered, Contiguous, objectSizes, "", totalSize );
    };

    /**
     * @brief Allocation
     * @param size number of objets in collection
     * @param totalSize total size of the collection
     * @param access type of access
     * @param objectSizes size of objects (constant or variable)
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value && MemoryType == Contiguous,
                             bool >::type
    allocate( const std::vector< std::vector< ValueType > > &values,
              JeveuxCollectionObjectSizes objectSizes = Variable ) {

        ASTERINTEGER totalSize = 0;
        for ( auto &val : values ) {
            totalSize += val.size();
        }

        genericAllocation( values.size(), Numbered, Contiguous, objectSizes, "", totalSize );

        for ( auto &val : values ) {
            push_back( val );
        }

        return true;
    };

    /**
     * @brief Allocation of one object by name
     */
    JeveuxCollObjValType allocateObject( const std::string &name, const ASTERINTEGER &nbValues ) {
        if ( size() > capacity() ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( size() ) + " vs " +
                      std::to_string( capacity() ) );
        }

        if ( !isNamed() ) {
            AS_ABORT( "Collection is not named" );
        }

        auto newIndex = getNewIndex();
        JeveuxCollObjValType obj( _name, newIndex, name, nbValues );

        _mapNumObject[std::string( strip( name.c_str() ) )] = newIndex;
        _indirection.push_back( newIndex );

        ++_objectCount;

        return obj;
    };

    /**
     * @brief Allocation of one object by name
     */
    JeveuxCollObjValType allocateObject( const std::string &name,
                                         const std::vector< ValueType > &values ) {
        auto obj = allocateObject( name, values.size() );
        obj->setValues( values );

        return obj;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    JeveuxCollObjValType allocateObject( const ASTERINTEGER &index, const ASTERINTEGER &nbValues ) {
        if ( contains( index ) ) {
            AS_ABORT( "Index already existing: " + std::to_string( index ) );
        }

        if ( size() > capacity() ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( size() ) + " vs " +
                      std::to_string( capacity() ) );
        }

        if ( !isNumeroted() ) {
            AS_ABORT( "Collection is not numeroted" );
        }

        JeveuxCollObjValType obj( _name, index, nbValues );
        _indirection.push_back( index );

        ++_objectCount;

        return obj;
    };

    /**
     * @brief Fast Allocation of one object at the end of a existing collection (contiguous)
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value && MemoryType == Contiguous,
                             void >::type
    fastAllocateObject( const ASTERINTEGER &index, const ASTERINTEGER &nbValues ) {
        if ( contains( index ) ) {
            AS_ABORT( "Index already existing: " + std::to_string( index ) );
        }

        if ( size() > capacity() ) {
            AS_ABORT( "Out of collection bounds: " + std::to_string( size() ) + " vs " +
                      std::to_string( capacity() ) );
        }

        if ( !isNumeroted() ) {
            AS_ABORT( "Collection is not numeroted" );
        }
        ASTERINTEGER taille = nbValues, ibid = index;

        std::string nameOfObject( "" );
        ValueType *valuePtr;

        CALLO_JUCROC_WRAP( _name, nameOfObject, &ibid, &taille, (void *)( &valuePtr ) );
        _indirection.push_back( index );

        ++_objectCount;
    };

    /**
     * @brief Fast access to start of a existing collection (contiguous)
     */
    template < typename T1 = AccessType, typename = IsSame< T1, ASTERINTEGER > >
    typename std::enable_if< std::is_same< T1, ASTERINTEGER >::value && MemoryType == Contiguous,
                             ValueType * >::type
    startPointer() {
        ValueType *valuePtr;
        std::string charJeveuxName( 32, ' ' );
        ASTERINTEGER num = 1;
        CALLO_JEXNUM( charJeveuxName, _name, &num );
        const std::string read( "L" );
        CALLO_JEVEUOC( charJeveuxName, read, (void *)( &valuePtr ) );

        return valuePtr;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    JeveuxCollObjValType allocateObject( const ASTERINTEGER &index,
                                         const std::vector< ValueType > &values ) {
        auto obj = allocateObject( index, values.size() );
        obj->setValues( values );

        return obj;
    };

    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    void push_back( const std::vector< ValueType > &values ) {
        auto obj = allocateObject( getNewIndex(), values.size() );
        obj->setValues( values );
    };
    /**
     * @brief Allocation of one object at the end of a existing collection
     */
    void push_back( const JeveuxCollObjValType &values ) {
        values->updateValuePointer();
        auto obj = allocateObject( getNewIndex(), values->size() );
        obj->setValues( *values );
    };
    /**
     * @brief Allocation of one object by name
     */
    void push_back( const std::string &name, const std::vector< ValueType > &values ) {
        auto obj = allocateObject( name, values.size() );
        obj->setValues( values );
    };

    /**
     * @brief Deallocate collection
     */
    void deallocate() {
        if ( _name != "" && get_sh_jeveux_status() == 1 )
            CALLO_JEDETR( _name );

        _capacity = 0;
        _objectCount = 0;
        _indirection.clear();
    };

    /**
     * @brief Methode permettant de construire une collection a partir d'une
     * collection existante en memoire Jeveux
     * @return Renvoit true si la construction s'est bien deroulee
     */
    bool build( bool force = false );

    /**
     * @brief Methode verifiant l'existence d'un objet de collection dans la
     * collection
     * @param name Chaine contenant le nom de l'objet
     * @return Renvoit true si l'objet existe dans la collection
     */
    bool contains( const std::string &name ) const;

    /**
     * @brief Methode verifiant l'existence d'un objet de collection dans la
     * collection
     * @param number entier
     * @return Renvoit true si l'objet existe dans la collection
     */
    bool contains( const ASTERINTEGER &number ) const;

    void updateValuePointer() const {
        const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )->build(
            true );
    }

    /**
     * @brief Methode permettant d'obtenir la liste des objets nommés dans la
     * collection
     * @return vecteur de noms d'objet de collection
     */
    std::vector< JeveuxChar32 > getObjectsNames() const;

    VectorLong getObjectsIndicies() const {
        VectorLong ret;
        ret.reserve( size() );

        for ( int key = 1; key <= size(); ++key ) {
            ret.push_back( key );
        }

        return ret;
    };

    ASTERINTEGER getObjectSize( const ASTERINTEGER &num ) const {
        ASTERINTEGER num2 = num;
        ASTERINTEGER valTmp;
        JeveuxChar8 param( "LONMAX" );
        std::string charval = std::string( 32, ' ' );
        std::string charJeveuxName( 32, ' ' );
        CALLO_JEXNUM( charJeveuxName, _name, &num2 );
        CALLO_JELIRA( charJeveuxName, param, &valTmp, charval );
        return valTmp;
    };

    std::vector< JeveuxCollObjValType > getObjects() const {
        std::vector< JeveuxCollObjValType > ret;
        ret.reserve( size() );

        for ( int key = 1; key <= size(); ++key ) {
            ret.push_back( this->operator[]( key ) );
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryCollectionObject = JeveuxCollObjValType();
        }

        return ret;
    };

    inline const JeveuxCollObjValType &operator[]( const ASTERINTEGER &position ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_isBuilt ) {
            AS_ABORT( "Collection not built: " + _name );
        }

#endif
        const auto lastPos = _indirection[_objectCount - 1];
        if ( position >= lastPos + 1 ) {
            auto obj = JeveuxCollObjValType();
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryEndCollectionObject = obj;
            return _temporaryEndCollectionObject;
        }
        if ( _temporaryCollectionObject.operator->() == nullptr ) {
            auto obj = JeveuxCollObjValType( _name, position, _isNamed );
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryCollectionObject = obj;
        } else {
            if ( _temporaryCollectionObject.operator->().use_count() > 1 ) {
                throw std::runtime_error( "Collection: " + _name + " index: " +
                                          std::to_string( position ) + " too many access" );
            } else {
                auto obj = JeveuxCollObjValType( _name, position, _isNamed );
                const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                    ->_temporaryCollectionObject = obj;
            }
        }
        _temporaryCollectionObject->updateValuePointer();
        return _temporaryCollectionObject;
    };

    inline JeveuxCollObjValType &operator[]( const ASTERINTEGER &position ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_isBuilt ) {
            AS_ABORT( "Collection not built: " + _name );
        }

#endif
        const auto lastPos = _indirection[_objectCount - 1];
        if ( position >= lastPos + 1 ) {
            auto obj = JeveuxCollObjValType();
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryEndCollectionObject = obj;
            return _temporaryEndCollectionObject;
        }
        if ( _temporaryCollectionObject.operator->() == nullptr ) {
            auto obj = JeveuxCollObjValType( _name, position, _isNamed );
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryCollectionObject = obj;
        } else {
            if ( _temporaryCollectionObject.operator->().use_count() > 1 ) {
                throw std::runtime_error( "Collection: " + _name + " index: " +
                                          std::to_string( position ) + " too many access" );
            } else {
                auto obj = JeveuxCollObjValType( _name, position, _isNamed );
                const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                    ->_temporaryCollectionObject = obj;
            }
        }
        _temporaryCollectionObject->updateValuePointer();
        return _temporaryCollectionObject;
    };

    inline const JeveuxCollObjValType &operator[]( const std::string &name ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_isBuilt ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( !_isNamed ) {
            AS_ABORT( "Collection " + _name + " is not named" );
        }
#endif
        const auto &curIter = _mapNumObject.find( strip( name ) );
        if ( curIter == _mapNumObject.end() ) {
            AS_ABORT( "Name not in collection: " + name );
        }
        const auto position = curIter->second;
        const auto lastPos = _indirection[_objectCount - 1];
        if ( position >= lastPos + 1 ) {
            auto obj = JeveuxCollObjValType();
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryEndCollectionObject = obj;
            return _temporaryEndCollectionObject;
        }
        if ( _temporaryCollectionObject.operator->() == nullptr ) {
            auto obj = JeveuxCollObjValType( _name, position, _isNamed );
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryCollectionObject = obj;
        } else {
            if ( _temporaryCollectionObject.operator->().use_count() > 1 ) {
                throw std::runtime_error( "Collection: " + _name + " index: " +
                                          std::to_string( position ) + " too many access" );
            } else {
                auto obj = JeveuxCollObjValType( _name, position, _isNamed );
                const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                    ->_temporaryCollectionObject = obj;
            }
        }
        _temporaryCollectionObject->updateValuePointer();
        return _temporaryCollectionObject;
    };

    inline JeveuxCollObjValType &operator[]( const std::string &name ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_isBuilt ) {
            AS_ABORT( "Collection not built: " + _name );
        }

        if ( !_isNamed ) {
            AS_ABORT( "Collection " + _name + " is not named" );
        }
#endif
        const auto &curIter = _mapNumObject.find( strip( name ) );
        if ( curIter == _mapNumObject.end() ) {
            AS_ABORT( "Name not in collection: " + name );
        }
        const auto position = curIter->second;
        const auto lastPos = _indirection[_objectCount - 1];
        if ( position >= lastPos + 1 ) {
            auto obj = JeveuxCollObjValType();
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryEndCollectionObject = obj;
            return _temporaryEndCollectionObject;
        }
        if ( _temporaryCollectionObject.operator->() == nullptr ) {
            auto obj = JeveuxCollObjValType( _name, position, _isNamed );
            const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                ->_temporaryCollectionObject = obj;
        } else {
            if ( _temporaryCollectionObject.operator->().use_count() > 1 ) {
                throw std::runtime_error( "Collection: " + _name + " index: " +
                                          std::to_string( position ) + " too many access" );
            } else {
                auto obj = JeveuxCollObjValType( _name, position, _isNamed );
                const_cast< JeveuxCollectionClass< ValueType, AccessType, MemoryType > * >( this )
                    ->_temporaryCollectionObject = obj;
            }
        }
        _temporaryCollectionObject->updateValuePointer();
        return _temporaryCollectionObject;
    };

    inline FastCollectionObject fastAccess( const int &pos ) const {
        return FastCollectionObject( *this, pos );
    };

    inline FastCollectionObject fastAccess( const std::string &name ) const {
        const auto &curIter = _mapNumObject.find( strip( name ) );
        if ( curIter == _mapNumObject.end() ) {
            AS_ABORT( "Name not in collection: " + name );
        }
        const auto pos = curIter->second;
        return fastAccess( pos );
    };

    inline ASTERINTEGER size() const { return _objectCount; };

    inline ASTERINTEGER capacity() const { return _capacity; };

    inline ASTERINTEGER totalSize() const {
        if ( exists() ) {
            ASTERINTEGER size;
            JeveuxChar8 param( "LONT" );
            std::string charval( 32, ' ' );
            CALLO_JELIRA( _name, param, &size, charval );
            return size;
        }
        return 0;
    };

    inline ASTERBOOL isNamed() const { return _isNamed; };

    inline ASTERBOOL isNumeroted() const { return ( !_isNamed ); };

    inline bool isContiguous() const {
        ASTERINTEGER nothin;
        JeveuxChar8 param( "STOCKAGE" );
        std::string charval( 32, ' ' );
        CALLO_JELIRA( _name, param, &nothin, charval );
        if ( strip( charval ) == "CONTIG" ) {
            return true;
        } else {
            return false;
        }
    };

    inline bool isBuilt() const { return _isBuilt; };

    /**
     * @brief Surcharge de l'operateur =
     */
    JeveuxCollectionClass &
    operator=( JeveuxCollectionClass< ValueType, AccessType, MemoryType > &toCopy ) {

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

        for ( int i = 1; i <= _objectCount; ++i ) {
            const auto curPos = _indirection[i - 1];
            const auto obj1 = this->operator[]( curPos );
            if ( !toCompar.contains( i ) ) {
                return false;
            }
            const auto obj2 = this->operator[]( i );
            if ( *obj1 != obj2 ) {
                return false;
            }
        }

        return true;
    };

    JeveuxCollectionClass< ValueType, AccessType > &operator*=( const ValueType &scal ) {
        this->build( false );
        auto size = this->size();

        for ( int i = 1; i <= _objectCount; ++i ) {
            const auto curPos = _indirection[i - 1];
            const auto obj1 = this->operator[]( curPos );
            ( *obj1 ) *= scal;
        }

        return *this;
    };

    /** @brief overload << operator */
    friend std::ostream &
    operator<<( std::ostream &os,
                const JeveuxCollectionClass< ValueType, AccessType, MemoryType > &toPrint ) {
        os << "JeveuxCollection: " << toPrint.getName() << "\n";
        os << "Size: " << std::to_string( toPrint.size() )
           << ", and capacity: " << std::to_string( toPrint.capacity() ) << ".\n";

        const auto size = toPrint.size();
        os << "List of objects: \n";
        for ( auto &[key, obj] : toPrint ) {
            obj->updateValuePointer();
            os << *obj << " \n";
        }

        return os;
    }
};

template < class ValueType, class AccessType, JeveuxCollectionMemoryStorageType MemoryType >
bool JeveuxCollectionClass< ValueType, AccessType, MemoryType >::build( bool force ) {

    if ( !exists() ) {
        return false;
    }

    ASTERINTEGER nbColObj, valTmp;
    JeveuxChar8 param( "NUTIOC" );
    std::string charval( 32, ' ' );
    CALLO_JELIRA( _name, param, &nbColObj, charval );

#ifdef ASTER_DEBUG_CXX
    // This is very stange that the size is 0
    // The most probable is that NUTIOC has not been setted
    if ( nbColObj <= 0 ) {
#ifdef ASTER_HAVE_MPI
        std::cout << getName() + " seems empty" << std::endl;
#else
        AS_ABORT( getName() + " seems empty" );
#endif
    }
#endif

    JeveuxChar8 param2( "NMAXOC" );
    CALLO_JELIRA( _name, param2, &_capacity, charval );

    if ( !force && _isBuilt && size() == nbColObj ) {
        return true;
    }

    param = "ACCES ";
    charval = std::string( 32, ' ' );
    CALLO_JELIRA( _name, param, &valTmp, charval );
    const std::string resu( charval, 0, 2 );

    if ( resu == "NO" ) {
        _isNamed = true;
    }

    _objectCount = 0;
    _indirection.clear();
    for ( ASTERINTEGER i = 1; i <= _capacity; ++i ) {
        if ( contains( i ) ) {
            _indirection.push_back( i );
            ++_objectCount;

            if ( _isNamed ) {
                std::string charJeveuxName( 32, ' ' );
                ASTERINTEGER num = i;
                CALLO_JEXNUM( charJeveuxName, _name, &num );
                const std::string objName = strip( charJeveuxName );
                _mapNumObject[objName] = i;
            }
        }
    }

    _isBuilt = true;
    return true;
};

template < class ValueType, class AccessType, JeveuxCollectionMemoryStorageType MemoryType >
bool JeveuxCollectionClass< ValueType, AccessType, MemoryType >::contains(
    const std::string &name ) const {
    std::string charJeveuxName( 32, ' ' );
    ASTERINTEGER returnBool;
    CALLO_JEXNOM( charJeveuxName, _name, name );
    CALLO_JEEXIN( charJeveuxName, &returnBool );
    if ( returnBool == 0 )
        return false;
    return true;
};

template < class ValueType, class AccessType, JeveuxCollectionMemoryStorageType MemoryType >
bool JeveuxCollectionClass< ValueType, AccessType, MemoryType >::contains(
    const ASTERINTEGER &number ) const {
    std::string charJeveuxName( 32, ' ' );
    ASTERINTEGER returnBool = number;
    CALLO_JEXNUM( charJeveuxName, _name, &returnBool );
    CALLO_JEEXIN( charJeveuxName, &returnBool );

    if ( returnBool == 0 )
        return false;
    return true;
};

template < class ValueType, class AccessType, JeveuxCollectionMemoryStorageType MemoryType >
std::vector< JeveuxChar32 >
JeveuxCollectionClass< ValueType, AccessType, MemoryType >::getObjectsNames() const {
    std::vector< JeveuxChar32 > toReturn;
    toReturn.reserve( size() );
    for ( ASTERINTEGER i = 1; i <= _objectCount; ++i ) {
        const auto obj = this->operator[]( i );
        toReturn.push_back( obj->getName() );
    }
    return toReturn;
};

/**
 * @class JeveuxCollection
 * @brief Enveloppe d'un pointeur intelligent vers un JeveuxCollectionClass
 * @author Nicolas Sellenet
 */
template < class ValueType, class AccessType = ASTERINTEGER,
           JeveuxCollectionMemoryStorageType MemoryType = Sparse >
class JeveuxCollection {
  public:
    typedef std::shared_ptr< JeveuxCollectionClass< ValueType, AccessType, MemoryType > >
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
              std::make_shared< JeveuxCollectionClass< ValueType, AccessType, MemoryType > >(
                  nom ) ) {}

    /**
     * @brief Constructeur dans le cas où AccessType est un NamesMap
     * @param name Chaine representant le nom de la collection
     */
    template < typename T1 = AccessType, typename = IsNotSame< T1, ASTERINTEGER > >
    JeveuxCollection( const std::string &nom, AccessType ptr )
        : _jeveuxCollectionPtr(
              std::make_shared< JeveuxCollectionClass< ValueType, AccessType, MemoryType > >(
                  nom, ptr ) ) {};

    ~JeveuxCollection() {};

    JeveuxCollection &
    operator=( const JeveuxCollection< ValueType, AccessType, MemoryType > &tmp ) {
        _jeveuxCollectionPtr = tmp._jeveuxCollectionPtr;
        return *this;
    };

    const JeveuxCollectionTypePtr &operator->() const { return _jeveuxCollectionPtr; };

    JeveuxCollectionClass< ValueType, AccessType, MemoryType > &operator*( void ) const {
        return *_jeveuxCollectionPtr;
    };

    bool exists() const {
        if ( _jeveuxCollectionPtr == nullptr )
            return false;
        if ( _jeveuxCollectionPtr.use_count() == 0 )
            return false;

        return _jeveuxCollectionPtr->exists();
    }

    auto begin() const { return _jeveuxCollectionPtr->begin(); };

    auto end() const { return _jeveuxCollectionPtr->end(); };

    auto cbegin() const { return _jeveuxCollectionPtr->cbegin(); };

    auto cend() const { return _jeveuxCollectionPtr->cend(); };
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

/** @typedef Definition d'une collection de type entier long */
typedef JeveuxCollection< ASTERINTEGER, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionLong;
/** @typedef Definition d'une collection de type entier court */
typedef JeveuxCollection< ASTERINTEGER4, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionShort;
/** @typedef Definition d'une collection de type double */
typedef JeveuxCollection< ASTERDOUBLE, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionReal;
/** @typedef Definition d'une collection de type double complex */
typedef JeveuxCollection< ASTERCOMPLEX, ASTERINTEGER, Contiguous >
    JeveuxContiguousCollectionComplex;
/** @typedef Definition d'une collection de JeveuxChar8 */
typedef JeveuxCollection< JeveuxChar8, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionChar8;
/** @typedef Definition d'une collection de JeveuxChar16 */
typedef JeveuxCollection< JeveuxChar16, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionChar16;
/** @typedef Definition d'une collection de JeveuxChar24 */
typedef JeveuxCollection< JeveuxChar24, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionChar24;
/** @typedef Definition d'une collection de JeveuxChar32 */
typedef JeveuxCollection< JeveuxChar32, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionChar32;
/** @typedef Definition d'une collection de JeveuxChar80 */
typedef JeveuxCollection< JeveuxChar80, ASTERINTEGER, Contiguous > JeveuxContiguousCollectionChar80;

#endif /* JEVEUXCOLLECTION_H_ */
