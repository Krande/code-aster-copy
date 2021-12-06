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
#include "Utilities/Tools.h"

ContactZone::ContactZone( const std::string name, const ModelPtr model )
    : DataStructure( name, 8, "CHAR_CONT_ZONE" ), _model( model ), _verbosity( 1 ),
      _checkNormal( true ) {
    // model has to be mechanics
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

bool ContactZone::build() {
    const auto mesh = getMesh();

    // check that there is no nodes in common
    auto slaveNodes = mesh->getNodesFromCells( getSlaveGroupOfCells() );
    auto masterNodes = mesh->getNodesFromCells( getMasterGroupOfCells() );

    std::sort( slaveNodes.begin(), slaveNodes.end() );
    std::sort( masterNodes.begin(), masterNodes.end() );

    VectorLong commonNodes;

    std::set_intersection( slaveNodes.begin(), slaveNodes.end(), masterNodes.begin(),
                           masterNodes.end(), std::back_inserter( commonNodes ) );

    if ( commonNodes.size() > 0 ) {
        UTMESS( "F", "CONTACT1_1" );
    }

    // check mesh orientation (normals)
    if ( checkNormals() ) {
        std::string slave = ljust( _slave, 24, ' ' );
        std::string master = ljust( _master, 24, ' ' );
        CALL_CHECKNORMALS( _model->getName().c_str(), slave.c_str(), master.c_str() );
    }
}
