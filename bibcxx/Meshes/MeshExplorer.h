#ifndef MESHEXPLORER_H_
#define MESHEXPLORER_H_

/**
 * @file MeshExplorer.h
 * @brief Fichier entete de la classe MeshExplorer
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

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"

class CellObject {
    const ASTERINTEGER _id;
    const VectorLong _listOfNodes;
    const ASTERINTEGER _type;

  public:
    CellObject( const ASTERINTEGER &id, const VectorLong &listOfNodes, const ASTERINTEGER &type )
        : _id( id ), _listOfNodes( listOfNodes ), _type( type ) {};

    ASTERINTEGER getNumberOfNodes() const { return _listOfNodes.size(); };

    const VectorLong &getNodes() const { return _listOfNodes; };

    ASTERINTEGER getId() const { return _id; };

    ASTERINTEGER getType() const { return _type; };

    /**
     * @brief
     */
    auto begin() const { return _listOfNodes.begin(); };

    /**
     * @brief
     */
    auto end() const { return _listOfNodes.end(); };

    /**
     * @brief
     */
    auto cbegin() const { return _listOfNodes.cbegin(); };

    /**
     * @brief
     */
    auto cend() const { return _listOfNodes.cend(); };
};

class CellsIteratorFromConnectivity {
  private:
    const JeveuxContiguousCollectionLong _connect;
    const JeveuxVectorLong _type;

  public:
    CellsIteratorFromConnectivity( const JeveuxContiguousCollectionLong &connect,
                                   const JeveuxVectorLong &type )
        : _connect( connect ), _type( type ) {
        _connect->build();
    };

    CellObject getCellObject( const int &pos ) const

    {
        const int size2 = _connect->size();
        if ( size2 <= 0 ) {
            AS_ABORT( "Connectivity not available" );
        }

        if ( pos > size2 || pos < 0 )
            return CellObject( -1, VectorLong(), -1 );
        auto &obj = ( *_connect )[pos + 1];
        obj->updateValuePointer();
        const ASTERINTEGER type = ( *_type )[pos];
        auto listNodes = obj->toVector();
        std::for_each( listNodes.begin(), listNodes.end(), []( ASTERINTEGER &d ) { d -= 1; } );
        return CellObject( pos, listNodes, type );
    };

    int size() const { return _type->size(); };
};

class CellsIteratorFromFiniteElementDescriptor {
  private:
    const JeveuxContiguousCollectionLong _connectAndType;

  public:
    CellsIteratorFromFiniteElementDescriptor( const JeveuxContiguousCollectionLong &connect )
        : _connectAndType( connect ) {
        _connectAndType->build();
    };

    CellObject getCellObject( const int &pos ) const

    {
        const int size2 = _connectAndType->size();
        if ( size2 <= 0 ) {
            AS_ABORT( "Connectivity not available" );
        }

        if ( pos > size2 || pos < 0 )
            return CellObject( -1, VectorLong(), -1 );
        auto &obj = ( *_connectAndType )[pos + 1];
        obj->updateValuePointer();
        const auto size = obj->size() - 1;
        const ASTERINTEGER type = ( *obj )[size];
        auto listNodes = obj->toVector();
        listNodes.pop_back();
        // Not yet zero-based.
        // std::for_each( listNodes.begin(), listNodes.end(), []( ASTERINTEGER &d ) { d -= 1; } );
        return CellObject( pos, listNodes, type );
    };

    int size() const { return _connectAndType->size(); };
};

/**
 * @class MeshExplorer
 * @brief Utility to loop over mesh cells
 * @author Nicolas Sellenet
 */
template < class ElemBuilder, typename... Args >
class MeshExplorer {
  private:
    const ElemBuilder _builder;

  public:
    /**
     * @brief Constructeur
     */
    MeshExplorer( const Args... a ) : _builder( a... ) {};

    /**
     * @brief Destructeur
     */
    ~MeshExplorer() {};

    struct const_iterator {
        int position;
        const ElemBuilder &builder;

        inline const_iterator( int memoryPosition, const ElemBuilder &test )
            : position( memoryPosition ), builder( test ) {};

        inline const_iterator( const const_iterator &iter )
            : position( iter.position ), builder( iter.builder ) {};

        inline const_iterator &operator=( const const_iterator &testIter ) {
            position = testIter.position;
            builder = testIter.builder;
            return *this;
        };

        inline const_iterator &operator++() {
            ++position;
            return *this;
        };

        inline bool operator==( const const_iterator &testIter ) const {
            return ( testIter.position == position );
        };

        inline bool operator!=( const const_iterator &testIter ) const {
            return ( testIter.position != position );
        };

        inline CellObject operator->() const { return builder.getCellObject( position ); };

        inline CellObject operator*() const { return builder.getCellObject( position ); };
    };

    inline CellObject operator[]( int position ) const {
        return _builder.getCellObject( position );
    };

    /**
     * @brief
     */
    const_iterator begin() const { return const_iterator( 0, _builder ); };

    /**
     * @brief
     * @todo revoir le fonctionnement du end car il peut provoquer de segfault
     */
    const_iterator end() const { return const_iterator( _builder.size(), _builder ); };

    /**
     * @brief Size of the explorer
     */
    ASTERINTEGER size() const { return _builder.size(); };
};

#endif /* MESHEXPLORER_H_ */
