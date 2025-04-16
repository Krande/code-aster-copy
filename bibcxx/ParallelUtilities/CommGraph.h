#ifndef COMMGRAPH_H_
#define COMMGRAPH_H_

/**
 * @file CommGraph.h
 * @brief Header of connection graph
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

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"
#include "ParallelUtilities/AsterMPI.h"

/**
 * @class CommGraph
 * @brief Class describing a communication graph
 * @author Nicolas Sellenet
 */
class CommGraph {
    /** @brief Vector of communication for the current process */
    VectorInt _commGraph;
    /** @brief Matchings find in graph */
    VectorPairInt _matchings;

  public:
    /** @brief Add an edge in comm graph */
    void addCommunication( const int &toAdd ) {
        const auto rank = getMPIRank();
        if ( rank == toAdd )
            throw std::runtime_error( "No self communication allowed" );
        _commGraph.push_back( toAdd );
    }

    /** @brief Get matchings */
    const VectorPairInt &getMatchings() const { return _matchings; };

    /** @brief Synchronize graph over processes and build matchings */
    void synchronizeOverProcesses();

    struct iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = PairInt;
        using pointer = PairInt *;
        using reference = PairInt &;

        iterator( pointer ptr ) : m_ptr( ptr ) {}

        reference operator*() const { return *m_ptr; }

        pointer operator->() const { return m_ptr; }

        iterator &operator++() {
            m_ptr++;
            // while ( *m_ptr == -1 )
            //     m_ptr++;
            return *this;
        }

        friend bool operator==( const iterator &a, const iterator &b ) {
            return a.m_ptr == b.m_ptr;
        };

        friend bool operator!=( const iterator &a, const iterator &b ) {
            return a.m_ptr != b.m_ptr;
        };

      private:
        pointer m_ptr;
    };

    iterator begin() { return iterator( &( _matchings[0] ) ); }

    iterator end() { return iterator( &( _matchings[_matchings.size()] ) ); }

    struct const_iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = PairInt;
        using pointer = PairInt *;
        using reference = PairInt &;

        const_iterator( pointer ptr ) : m_ptr( ptr ) {}

        reference operator*() const { return *m_ptr; }

        pointer operator->() const { return m_ptr; }

        const_iterator &operator++() {
            m_ptr++;
            while ( m_ptr->first == -1 )
                m_ptr++;
            return *this;
        }

        friend bool operator==( const const_iterator &a, const const_iterator &b ) {
            return a.m_ptr == b.m_ptr;
        };

        friend bool operator!=( const const_iterator &a, const const_iterator &b ) {
            return a.m_ptr != b.m_ptr;
        };

      private:
        pointer m_ptr;
    };

    const_iterator cbegin() { return const_iterator( &( _matchings[0] ) ); }

    const_iterator cend() { return const_iterator( &( _matchings[_matchings.size()] ) ); }
};

using CommGraphPtr = std::shared_ptr< CommGraph >;

#endif /* ASTER_HAVE_MPI */

#endif /* COMMGRAPH_H_ */
