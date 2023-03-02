/**
 * @file ParallelDOFNumbering.cxx
 * @brief Implementation de ParallelDOFNumbering
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Meshes/Joints.h"

#ifdef ASTER_HAVE_MPI

Joints::Joints() : Joints( DataStructureNaming::getNewName() ){};

Joints::Joints( const std::string name )
    : DataStructure( name, 19, "DOMJOINTS" ),
      _domj( JeveuxVectorLong( getName() + ".DOMJ" ) ),
      _send( JeveuxCollectionLong( getName() + ".SEND" ) ),
      _recv( JeveuxCollectionLong( getName() + ".RECV" ) ){};

void Joints::setOppositeDomains( const VectorLong &oppositeDomains ) {
    ( *_domj ) = oppositeDomains;
};

JeveuxVectorLong Joints::getOppositeDomains() const {
    if ( _domj->exists() ) {
        _domj->updateValuePointer();
    }

    return _domj;
}

bool Joints::build() {
    _send->build();
    _recv->build();

    return true;
}

void Joints::setSendedElements( const VectorOfVectorsLong &send ) {
    if ( send.size() > 0 ) {
        _send->allocateSparseNumbered( send.size() );

        ASTERINTEGER i = 1;
        for ( auto &send_i : send ) {
            if ( send_i.size() > 0 ) {
                auto obj = _send->allocateObject( i, send_i );
            }
            i++;
        }
    }
};

void Joints::setReceivedElements( const VectorOfVectorsLong &recv ) {
    if ( recv.size() > 0 ) {
        _recv->allocateSparseNumbered( recv.size() );

        ASTERINTEGER i = 1;
        for ( auto &recv_i : recv ) {
            if ( recv_i.size() > 0 ) {
                auto obj = _recv->allocateObject( i, recv_i );
            }
            i++;
        }
    }
};

#endif /* ASTER_HAVE_MPI */
