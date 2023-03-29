/**
 * @file IncompleteMesh.cxx
 * @brief Implementation de IncompleteMesh
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
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "Meshes/IncompleteMesh.h"

const std::map< int, std::set< int > > &IncompleteMesh::buildReverseConnectivity() {
    if ( _bReverseConnex )
        return _reverseConnex;
    const auto connex = getConnectivityExplorer();
    int elemId = 0;
    for ( const auto &element : connex ) {
        for ( const auto &nodeId : element ) {
            _reverseConnex[nodeId - 1].insert( elemId );
        }
        ++elemId;
    }
    _bReverseConnex = true;
    return _reverseConnex;
};

void IncompleteMesh::deleteReverseConnectivity() {
    // free memory
    _reverseConnex = std::map< int, std::set< int > >();
    _bReverseConnex = false;
};

#endif /* ASTER_HAVE_MPI */
