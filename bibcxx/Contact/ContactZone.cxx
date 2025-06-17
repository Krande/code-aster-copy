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
    : DSWithCppPickling( name, 8, "CHAR_CONT_ZONE" ),
      _model( nullptr ),
      _verbosity( 1 ),
      _checkNormal( true ),
      _smoothing( false ),
      _contParam( std::make_shared< ContactParameter >() ),
      _fricParam( std::make_shared< FrictionParameter >() ),
      _pairParam( std::make_shared< PairingParameter >() ),
      _meshPairing( std::make_shared< MeshPairing >( getName() + ".APMA" ) ) {};

ContactZone::ContactZone( const py::tuple &tup ) : ContactZone( tup[0].cast< std::string >() ) {
    int i = 0;
    _model = tup[++i].cast< ModelPtr >();
    _verbosity = tup[++i].cast< ASTERINTEGER >();
    _checkNormal = tup[++i].cast< bool >();
    _smoothing = tup[++i].cast< bool >();
    _contParam = tup[++i].cast< ContactParameterPtr >();
    _fricParam = tup[++i].cast< FrictionParameterPtr >();
    _pairParam = tup[++i].cast< PairingParameterPtr >();
    _meshPairing = tup[++i].cast< MeshPairingPtr >();
}
py::tuple ContactZone::_getState() const {
    return py::make_tuple( getName(), _model, _verbosity, _checkNormal, _smoothing, _contParam,
                           _fricParam, _pairParam, _meshPairing );
}

void ContactZone::setVerbosity( const ASTERINTEGER &level ) {
    _verbosity = level;
    _meshPairing->setVerbosity( getVerbosity() );
}

bool ContactZone::pairing( ASTERDOUBLE &dist_pairing, ASTERDOUBLE &pair_tole ) {
    return _meshPairing->compute( dist_pairing, pair_tole );
}

bool ContactZone::build( const ModelPtr model ) {
    _model = model;
    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );

    _meshPairing->setMesh( _model->getMesh() );
    _meshPairing->build();
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
