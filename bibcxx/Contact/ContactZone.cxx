/**
 * @file ContactZone.cxx
 * @brief Implementation de Contact
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "aster_fort_mesh.h"

#include "Contact/ContactZone.h"
#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

ContactZone::ContactZone( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONT_ZONE" ), _model( model ), _verbosity( 1 ),
      _checkNormal( true ), _contParam( boost::make_shared< ContactParameter >() ),
      _fricParam( boost::make_shared< FrictionParameter >() ),
      _pairParam( boost::make_shared< PairingParameter >() ) {
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

bool ContactZone::build() {
    auto mesh = getMesh();

    // check that there is no nodes in common
    auto slaveNodes_lc = mesh->getNodesFromCells( getSlaveGroupOfCells(), false, true );
    auto masterNodes_lc = mesh->getNodesFromCells( getMasterGroupOfCells(), false, true );

    VectorLong commonNodes;

    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        VectorLong slaveNodes_gl;
        AsterMPI::all_gather( slaveNodes_lc, slaveNodes_gl );
        commonNodes = set_intersection( slaveNodes_gl, masterNodes_lc );
#endif
    } else {
        commonNodes = set_intersection( slaveNodes_lc, masterNodes_lc );
    }

    ASTERINTEGER size_inter_gl = commonNodes.size();

    // share error
#ifdef ASTER_HAVE_MPI
    if ( mesh->isParallel() ) {
        ASTERINTEGER size_inter_lc = size_inter_gl;
        size_inter_gl = 0;
        AsterMPI::all_reduce( size_inter_lc, size_inter_gl, MPI_SUM );
    }
#endif

    if ( size_inter_gl > 0 ) {
        UTMESS( "F", "CONTACT1_1" );
    }

    // check mesh orientation (normals)
    if ( checkNormals() ) {
        std::string slave = ljust( getSlaveGroupOfCells(), 24, ' ' );
        std::string master = ljust( getMasterGroupOfCells(), 24, ' ' );
        CALL_CHECKNORMALS( _model->getName().c_str(), slave.c_str(), master.c_str() );
    }
    return true;
}
