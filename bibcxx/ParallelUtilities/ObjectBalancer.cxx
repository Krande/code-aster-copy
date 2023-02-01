/**
 * @file CommGraph.cxx
 * @brief Implementation of an object balancer
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

#include "ParallelUtilities/ObjectBalancer.h"

void ObjectBalancer::prepareCommunications() {
    _graph->synchronizeOverProcesses();
    const auto rank = getMPIRank();
    int tag = 0;
    for ( const auto proc : *_graph ) {
        VectorInt tmp( 1, -1 );
        if ( rank > proc ) {
            tmp[0] = _sendList[proc].size();
            AsterMPI::send( tmp, proc, tag );
            AsterMPI::receive( tmp, proc, tag );
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << rank << " received " << tmp[0] << " from #" << proc << std::endl;
#endif
            _recvSize[proc] = tmp[0];
        } else {
            AsterMPI::receive( tmp, proc, tag );
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << rank << " received " << tmp[0] << " from #" << proc << std::endl;
#endif
            _recvSize[proc] = tmp[0];
            tmp[0] = _sendList[proc].size();
            AsterMPI::send( tmp, proc, tag );
        }
        ++tag;
    }
};
