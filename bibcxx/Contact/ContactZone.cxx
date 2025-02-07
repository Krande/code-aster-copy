/**
 * @file ContactZone.cxx
 * @brief Implementation de Contact
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

#include "Contact/ContactZone.h"

#include "aster_fort_mesh.h"

#include "Meshes/MeshPairing.h"
#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

ContactZone::ContactZone( const std::string name )
    : DataStructure( name, 8, "CHAR_CONT_ZONE" ),
      _model( nullptr ),
      _verbosity( 1 ),
      _checkNormal( true ),
      _smoothing( false ),
      _contParam( std::make_shared< ContactParameter >() ),
      _fricParam( std::make_shared< FrictionParameter >() ),
      _pairParam( std::make_shared< PairingParameter >() ),
      _meshPairing( std::make_shared< MeshPairing >( getName() + ".APMA" ) ) {};

void ContactZone::setVerbosity( const ASTERINTEGER &level ) {
    _verbosity = level;
    if ( _meshPairing != nullptr ) {
        _meshPairing->setVerbosity( getVerbosity() );
    }
}

bool ContactZone::pairing( ASTERDOUBLE &dist_pairing, ASTERDOUBLE &pair_tole ) {
    return _meshPairing->compute( dist_pairing, pair_tole );
}

bool ContactZone::build( const ModelPtr model ) {
    _model = model;
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );

    _meshPairing->initObjects( _model->getMesh() );
    _meshPairing->setVerbosity( getVerbosity() );

    // Some checks
    auto hasCommonNodes = _meshPairing->hasCommonNodes();
    if ( hasCommonNodes ) {
        UTMESS( "F", "CONTACT1_1" );
    }

    if ( getContactParameter()->getAlgorithm() == ContactAlgo::Nitsche ) {
        auto surf2Volu = getSlaveCellsSurfToVolu();
        if ( surf2Volu.size() == 0 ) {
            UTMESS( "F", "CONTACT1_3" );
        }
    }

    // Check mesh orientation (normals)
    if ( checkNormals() ) {
        _meshPairing->checkNormals( _model );
    }

    return true;
}
