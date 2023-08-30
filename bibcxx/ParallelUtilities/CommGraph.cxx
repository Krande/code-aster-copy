/**
 * @file CommGraph.cxx
 * @brief Implementation de CommGraph
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

#include "ParallelUtilities/CommGraph.h"

#ifdef ASTER_HAVE_MPI

#include "ParallelUtilities/AsterMPI.h"

void CommGraph::synchronizeOverProcesses() {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorInt myComms( nbProcs, 0 );
    for ( const auto &curComm : _commGraph ) {
        myComms[curComm] = 1;
    }
    VectorInt outValues;
    AsterMPI::all_gather( myComms, outValues );
    for ( int iProc1 = 0; iProc1 < nbProcs; ++iProc1 ) {
        for ( int iProc2 = 0; iProc2 < nbProcs; ++iProc2 ) {
            const int posIt = iProc1 * nbProcs + iProc2;
            const int posIt2 = iProc2 * nbProcs + iProc1;
            if ( outValues[posIt] != 0 )
                outValues[posIt2] = 1;
            if ( outValues[posIt2] != 0 )
                outValues[posIt] = 1;
        }
    }
    int nbEdges = 0;
    for ( const auto &val : outValues ) {
        if ( val != 0 )
            ++nbEdges;
    }
    nbEdges /= 2;
    VectorInt mask( nbProcs * nbProcs, 0 );
    int nMatch = 1;
    while ( nbEdges > 0 ) {
        VectorInt marks( nbProcs, 0 );
        for ( int iProc1 = 0; iProc1 < nbProcs; ++iProc1 ) {
            for ( int iProc2 = 0; iProc2 < nbProcs; ++iProc2 ) {
                const int posIt = iProc1 * nbProcs + iProc2;
                if ( outValues[posIt] == 1 && marks[iProc1] == 0 && marks[iProc2] == 0 ) {
                    outValues[posIt] = 0;
                    mask[posIt] = nMatch;
                    const int posIt2 = iProc2 * nbProcs + iProc1;
                    outValues[posIt2] = 0;
                    mask[posIt2] = nMatch;
                    --nbEdges;
                    marks[iProc1] = 1;
                    marks[iProc2] = 1;
                }
            }
        }
        ++nMatch;
    }
    --nMatch;
    _matchings = VectorPairInt( nMatch, std::make_pair( -1, -1 ) );
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        auto num = mask[rank * nbProcs + iProc];
        if ( num > nMatch )
            throw std::runtime_error( "Error in graph building" );
        if ( num != 0 ) {
            _matchings[num - 1] = std::make_pair( num, iProc );
#ifdef ASTER_DEBUG_CXX
            std::cout << "In matching " << num << " #" << rank << " will communicate with #"
                      << iProc << std::endl;
#endif
        }
    }
};

#endif /* ASTER_HAVE_MPI */
