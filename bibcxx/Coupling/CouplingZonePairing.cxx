/**
 * @file CouplingZonePairing.cxx
 * @brief Implementation de Coupling
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "Coupling/CouplingZonePairing.h"

#include "aster_fort_ds.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Messages/Messages.h"
#include "Utilities/Tools.h"

CouplingZonePairing::CouplingZonePairing( const BaseMeshPtr mesh, const ASTERINTEGER verbosity )
    : DataStructure( ResultNaming::getNewResultName(), 8, "CHAR_CPL_ZONE" ),
      _meshPairing( std::make_shared< MeshPairing >( getName() + ".APMA" ) ),
      _algo( CouplingMethod::Undefined ),
      _coef_pena( 0. ) {

    _meshPairing->setMesh( mesh );
    _meshPairing->setVerbosity( verbosity );
    _meshPairing->setMethod( PairingMethod::Legacy );
};

void CouplingZonePairing::setSlaveGroupsOfCells( const VectorString &groupName ) {
    if ( groupName.size() != 1 ) {
        raiseAsterError( "Only one group allowed for the moment" );
    }
    _meshPairing->setSlaveGroupOfCells( groupName[0] );
};

void CouplingZonePairing::setMasterGroupsOfCells( const VectorString &groupName ) {
    if ( groupName.size() != 1 ) {
        raiseAsterError( "Only one group allowed for the moment" );
    }
    _meshPairing->setMasterGroupOfCells( groupName[0] );
};

void CouplingZonePairing::setVerbosity( const ASTERINTEGER verbosity ) {
    _meshPairing->setVerbosity( verbosity );
};

void CouplingZonePairing::check( const ModelPtr model ) const {
    // Some checks
    auto hasCommonNodes = _meshPairing->hasCommonNodes();
    if ( hasCommonNodes ) {
        UTMESS( "F", "CONTACT1_1" );
    }

    if ( _algo == CouplingMethod::Nitsche ) {
        auto surf2Volu = _meshPairing->getSlaveCellsSurfToVolu();
        if ( surf2Volu.size() == 0 ) {
            UTMESS( "F", "CONTACT1_3" );
        }
    }

    // Check mesh orientation (normals)
    _meshPairing->checkNormals( model );
}

bool CouplingZonePairing::build() { return _meshPairing->build(); };

void CouplingZonePairing::setPairingParameters( const PairingParameterPtr params ) {
    _params = params;

    _meshPairing->setMethod( _params->getPairingMethod() );
};

ASTERBOOL CouplingZonePairing::compute() {

    // Pairing
    const auto returnValue =
        _meshPairing->compute( _params->getDistanceRatio(), _params->getPairingTolerance(),
                               _params->getAreaIntersectionTolerance() );

    if ( !returnValue ) {
        return returnValue;
    }

    return true;
};

VectorOfVectorsReal
CouplingZonePairing::getIntersectionPoints( const CoordinatesSpace coorSpace ) const {

    VectorOfVectorsReal returnValue;
    const ASTERINTEGER nbPairs = this->getNumberOfPairs();
    if ( nbPairs == 0 ) {
        raiseAsterError( "No contact pairs: was the pairing performed correctly? " );
    }
    returnValue.reserve( nbPairs );
    for ( auto iPair = 0; iPair < nbPairs; iPair++ ) {
        const VectorOfVectorsReal interOnZone =
            _meshPairing->getIntersectionPoints( iPair, coorSpace );
        const ASTERINTEGER nbInter = _meshPairing->getNumberOfIntersectionPoints( iPair );
        if ( nbInter != 0 ) {
            VectorReal vectVale;
            for ( auto iInter = 0; iInter < nbInter; iInter++ ) {
                vectVale.insert( vectVale.end(), interOnZone[iInter].begin(),
                                 interOnZone[iInter].end() );
            };
            returnValue.push_back( vectVale );
        }
    }
    return returnValue;
};

VectorOfVectorsReal
CouplingZonePairing::getIntersectionPoints( const ASTERINTEGER indexPair,
                                            const CoordinatesSpace coorSpace ) const {

    return _meshPairing->getIntersectionPoints( indexPair, coorSpace );
};
