/**
 * @file MeshPairing.cxx
 * @brief Implementation of MeshPairing class
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

#include "Meshes/MeshPairing.h"

#include "aster_fort_ds.h"
#include "aster_fort_mesh.h"

#include "Messages/Messages.h"
#include "Utilities/Tools.h"

MeshPairing::MeshPairing( const std::string name )
    : DSWithCppPickling( name, 13, "MESH_PAIRING" ),
      _mesh( nullptr ),
      _masterInverseConnectivity( JeveuxCollectionLong( getName() + ".CM" ) ),
      _slaveInverseConnectivity( JeveuxCollectionLong( getName() + ".CS" ) ),
      _masterNeighbors( JeveuxCollectionLong( getName() + ".MN" ) ),
      _slaveNeighbors( JeveuxCollectionLong( getName() + ".SN" ) ),
      _verbosity( 1 ),
      _slaveCellsGroup( " " ),
      _masterCellsGroup( " " ),
      _zoneHaveBeenDefined( false ),
      _nbPairs( 0 ),
      _method( PairingMethod::Fast ),
      _currentCoordinates( nullptr ) {};

MeshPairing::MeshPairing( const py::tuple &tup ) : MeshPairing( tup[0].cast< std::string >() ) {
    int i = 0;
    setMesh( tup[++i].cast< BaseMeshPtr >() );
    _verbosity = tup[++i].cast< ASTERINTEGER >();
    _method = tup[++i].cast< PairingMethod >();
    std::string slavGrp = tup[++i].cast< std::string >();
    std::string masterGrp = tup[++i].cast< std::string >();
    _nbPairs = tup[++i].cast< ASTERINTEGER >();
    _pairs = tup[++i].cast< VectorLong >();
    _nbPoinInte = tup[++i].cast< VectorLong >();
    _poinInteSlav = tup[++i].cast< VectorReal >();
    _slavSurf2Volu = tup[++i].cast< MapLong >();
    setPair( slavGrp, masterGrp );
};

py::tuple MeshPairing::_getState() const {
    return py::make_tuple( getName(), _mesh, _verbosity, _method, _slaveCellsGroup,
                           _masterCellsGroup, _nbPairs, _pairs, _nbPoinInte, _poinInteSlav,
                           _slavSurf2Volu );
};

void MeshPairing::setMesh( const BaseMeshPtr &mesh ) {
    _mesh = mesh;
    if ( _currentCoordinates == nullptr ) {
        _currentCoordinates =
            std::make_shared< MeshCoordinatesField >( *( _mesh->getCoordinates() ) );
    }
}

void MeshPairing::setSlaveGroupOfCells( const std::string &groupName ) {
    _slaveCellsGroup = groupName;
};

void MeshPairing::setMasterGroupOfCells( const std::string &groupName ) {
    _masterCellsGroup = groupName;
};

void MeshPairing::setPair( const std::string &groupNameSlav, const std::string &groupNameMast ) {

    MeshPairing::setMasterGroupOfCells( groupNameMast );
    MeshPairing::setSlaveGroupOfCells( groupNameSlav );

    this->build();
};

ASTERBOOL
MeshPairing::buildInverseConnectivity() {
    // Create master inverse connectivity
    ASTERINTEGER nbMaster = getMasterCells().size();
    std::string base( "G" );

    VectorLong masterCells;
    masterCells.reserve( _masterCells.size() );

    // Shifting for fortran
    for ( auto cell : _masterCells )
        masterCells.push_back( cell + 1 );

    _masterInverseConnectivity->deallocate();
    CALL_CNCINV( getMesh()->getName().c_str(), masterCells.data(), &nbMaster, base.c_str(),
                 _masterInverseConnectivity->getName().c_str() );
    _masterInverseConnectivity->build();

    // create slave inverse connectivity
    ASTERINTEGER nbSlave = _slaveCells.size();

    VectorLong slaveCells;
    slaveCells.reserve( _slaveCells.size() );

    // Shifting for fortran
    for ( auto cell : _slaveCells )
        slaveCells.push_back( cell + 1 );

    _slaveInverseConnectivity->deallocate();
    CALL_CNCINV( getMesh()->getName().c_str(), slaveCells.data(), &nbSlave, base.c_str(),
                 _slaveInverseConnectivity->getName().c_str() );
    _slaveInverseConnectivity->build();

    return true;
}

ASTERBOOL MeshPairing::buildCellsNeighbors() {
    if ( _masterNeighbors->exists() && _slaveNeighbors->exists() ) {
        return true;
    }

    ASTERINTEGER ind_max, ind_min;

    // Get master neighbors
    ASTERINTEGER nbMaster = getMasterCells().size();
    if ( nbMaster > 0 ) {
        ind_max = *std::max_element( _masterCells.begin(), _masterCells.end() ) + 1;
        ind_min = *std::min_element( _masterCells.begin(), _masterCells.end() ) + 1;

        std::string invmcn_name = ljust( _masterInverseConnectivity->getName(), 24, ' ' );
        std::string mn_name = ljust( _masterNeighbors->getName(), 24, ' ' );

        VectorLong masterCells;
        masterCells.reserve( _masterCells.size() );
        for ( auto cell : _masterCells )
            masterCells.push_back( cell + 1 );

        _masterNeighbors->deallocate();
        CALL_CNVOIS( getMesh()->getName(), masterCells.data(), invmcn_name, &nbMaster, &ind_min,
                     &ind_max, mn_name );

        _masterNeighbors->build();
    }

    // Get slave neighbors
    ASTERINTEGER nbSlave = getSlaveCells().size();
    if ( nbSlave > 0 ) {
        ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() ) + 1;
        ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() ) + 1;

        std::string invscn_name = ljust( _slaveInverseConnectivity->getName(), 24, ' ' );
        std::string sn_name = ljust( _slaveNeighbors->getName(), 24, ' ' );

        VectorLong slaveCells;
        slaveCells.reserve( _slaveCells.size() );
        for ( auto cell : _slaveCells )
            slaveCells.push_back( cell + 1 );
        CALL_CNVOIS( getMesh()->getName(), slaveCells.data(), invscn_name, &nbSlave, &ind_min,
                     &ind_max, sn_name );

        _slaveNeighbors->build();
    }

    return true;
}

void MeshPairing::clearResult() {
    _nbPairs = 0;
    _pairs.clear();
    _nbPoinInte.clear();
    _poinInteSlav.clear();
}

ASTERBOOL MeshPairing::surfacesHasBeenDefined() {
    ASTERBOOL returnValue;
    returnValue = ( _slaveCellsGroup.size() != 0 && _masterCellsGroup.size() != 0 );
    return returnValue;
}

ASTERBOOL MeshPairing::compute( ASTERDOUBLE &dist_pairing, ASTERDOUBLE &pair_tole ) {

    build();
    CALL_JEMARQ();

    // Get and define some input parameters
    VectorLong masterCells = getMasterCells();
    VectorLong masterNodes = getMasterNodes();
    VectorLong slaveCells = getSlaveCells();
    ASTERINTEGER nbCellMaster = masterCells.size();
    ASTERINTEGER nbNodeMaster = masterNodes.size();
    ASTERINTEGER nbCellSlave = slaveCells.size();
    ASTERINTEGER verbosity = getVerbosity();
    std::string mastConnexInveName = getMasterInverseConnectivityName();
    std::string mastNeighName = getMasterNeighName();
    std::string slavNeighName = getSlaveNeighName();

    // Update the numbering for fortran
    std::for_each( masterCells.begin(), masterCells.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( masterNodes.begin(), masterNodes.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( slaveCells.begin(), slaveCells.end(), []( ASTERINTEGER &d ) { d += 1; } );

    // Set pairs numbers to 0
    ASTERINTEGER nb_pairs = 0;

    // Clear pairing results
    this->clearResult();

    // Method
    ASTERINTEGER method;
    if ( _method == PairingMethod::Fast ) {
        method = 1;
    } else if ( _method == PairingMethod::Legacy ) {
        method = 2;
    } else if ( _method == PairingMethod::BrutForce ) {
        method = 3;
    } else {
        raiseAsterError( "Unknown method of pairing " );
    }

    // Main routine for pairing
    CALLO_PAIRWRAP( &method, _mesh->getName(), _currentCoordinates->getName(), mastConnexInveName,
                    mastNeighName, slavNeighName, &pair_tole, &dist_pairing, &verbosity,
                    &nbCellMaster, masterCells.data(), &nbCellSlave, slaveCells.data(),
                    &nbNodeMaster, masterNodes.data(), &nb_pairs, getBasename() );

    // Get number of pairs
    _nbPairs = nb_pairs;

    // Output JEVEUX objects
    if ( _nbPairs != 0 ) {
        if ( !_jvPairs.exists() ) {
            _jvPairs = JeveuxVectorLong( getPairsName() );
        }
        if ( !_jvNbInterPoints.exists() ) {
            _jvNbInterPoints = JeveuxVectorLong( getNbInterName() );
        }
        if ( !_jvInterSlavePoints.exists() ) {
            _jvInterSlavePoints = JeveuxVectorReal( getCoorInterName() );
        }
        _jvPairs->updateValuePointer();
        _jvNbInterPoints->updateValuePointer();
        _jvInterSlavePoints->updateValuePointer();

        // Set output values
        _pairs = _jvPairs->toVector();
        std::transform( _pairs.begin(), _pairs.end(), _pairs.begin(),
                        []( ASTERINTEGER &indexCell ) -> ASTERINTEGER { return --indexCell; } );
        _nbPoinInte = _jvNbInterPoints->toVector();
        _poinInteSlav = _jvInterSlavePoints->toVector();
    }

    CALL_JEDEMA();

    return true;
}

VectorPairLong MeshPairing::getListOfPairs() const {

    VectorPairLong returnValue;
    ASTERINTEGER nbPairs = getNumberOfPairs();

    if ( nbPairs == 0 ) {
        throw std::runtime_error( "No pairs !" );
    }

    if ( _pairs.size() == 0 ) {
        throw std::runtime_error( "No pairs from Fortran!" );
    }

    returnValue.reserve( nbPairs );

    for ( auto iPair = 0; iPair < nbPairs; iPair++ ) {
        returnValue.push_back( std::make_pair( _pairs[2 * iPair], _pairs[2 * iPair + 1] ) );
    }

    return returnValue;
}

ASTERINTEGER MeshPairing::getNumberOfIntersectionPoints( const ASTERINTEGER &indexPair ) const {

    ASTERINTEGER nbPairs = getNumberOfPairs();

    if ( indexPair < 0 || indexPair >= nbPairs ) {
        throw std::out_of_range( "The pair index should be between 0  and " +
                                 std::to_string( nbPairs ) );
    }

    return _nbPoinInte[indexPair];
}

VectorLong MeshPairing::getNumberOfIntersectionPoints() const {

    ASTERINTEGER nbPairs = getNumberOfPairs();

    VectorLong returnValue;

    for ( auto indexPair = 0; indexPair < nbPairs; indexPair++ ) {
        ASTERINTEGER nbInter = getNumberOfIntersectionPoints( indexPair );
        returnValue.push_back( nbInter );
    }
    return returnValue;
}

VectorOfVectorsReal MeshPairing::getIntersectionPoints( const ASTERINTEGER &indexPair,
                                                        const CoordinatesSpace coorSpace ) const {

    ASTERINTEGER nbPairs = getNumberOfPairs();

    if ( indexPair < 0 || indexPair >= nbPairs ) {
        throw std::out_of_range( "The pair index should be between 0  and " +
                                 std::to_string( nbPairs ) );
    }

    ASTERINTEGER nbInter = getNumberOfIntersectionPoints( indexPair );

    VectorOfVectorsReal returnValue;

    if ( coorSpace == CoordinatesSpace::Global ) {
        VectorReal coorGlobPoint;
        coorGlobPoint.reserve( 3 );
        coorGlobPoint.push_back( 0.0 );
        coorGlobPoint.push_back( 0.0 );
        coorGlobPoint.push_back( 0.0 );

        // Objects for output
        VectorReal poinInteReal;
        poinInteReal.reserve( 3 * 8 );
        CALLO_INTEPOINCOORWRAP( _mesh->getName(), _currentCoordinates->getName(), getBasename(),
                                &indexPair, poinInteReal.data() );

        if ( nbInter == 0 ) {
            throw std::out_of_range( "No intersection in pair " + std::to_string( indexPair ) );
        } else {
            returnValue.reserve( 3 * nbInter );
            for ( auto iInter = 0; iInter < nbInter; iInter++ ) {
                coorGlobPoint[0] = poinInteReal[3 * iInter];
                coorGlobPoint[1] = poinInteReal[3 * iInter + 1];
                coorGlobPoint[2] = poinInteReal[3 * iInter + 2];
                returnValue.push_back( coorGlobPoint );
            }
        }
    } else if ( coorSpace == CoordinatesSpace::Slave ) {

        VectorReal coorParaPoint;
        coorParaPoint.reserve( 2 );
        coorParaPoint.push_back( 0.0 );
        coorParaPoint.push_back( 0.0 );

        if ( nbInter != 0 ) {
            returnValue.reserve( nbInter );
            for ( auto iInter = 0; iInter < nbInter; iInter++ ) {
                coorParaPoint[0] = _poinInteSlav[16 * ( indexPair ) + iInter];
                coorParaPoint[1] = _poinInteSlav[16 * ( indexPair ) + iInter + 8];
                returnValue.push_back( coorParaPoint );
            }
        }
    }
    return returnValue;
}

ASTERDOUBLE MeshPairing::getIntersectionArea( const ASTERINTEGER &indexPair ) const {

    ASTERINTEGER nbPairs = getNumberOfPairs();

    if ( indexPair < 0 || indexPair >= nbPairs ) {
        throw std::out_of_range( "The pair index should be between 0  and " +
                                 std::to_string( nbPairs ) );
    }

    double coorPoint;
    VectorReal poinInte;
    ASTERDOUBLE inteArea = 0.0;

    ASTERINTEGER nbInter = getNumberOfIntersectionPoints( indexPair );
    if ( nbInter != 0 ) {
        for ( auto iInter = 0; iInter < nbInter; iInter++ ) {
            coorPoint = _poinInteSlav[16 * indexPair + iInter];
            poinInte.push_back( coorPoint );
            coorPoint = _poinInteSlav[16 * indexPair + 8 + iInter];
            poinInte.push_back( coorPoint );
        }
        CALLO_INTECELLAREAWRAP( _mesh->getName(), &nbInter, poinInte.data(), &inteArea );
    }
    return inteArea;
}

std::vector< VectorReal > MeshPairing::getQuadraturePoints( const ASTERINTEGER &indexPair ) const {

    std::vector< VectorReal > returnValue;

    ASTERINTEGER nbPairs = getNumberOfPairs();
    if ( indexPair < 0 || indexPair >= nbPairs ) {
        throw std::out_of_range( "The pair index should be between 0  and " +
                                 std::to_string( nbPairs ) );
    }

    VectorReal coorPoint;
    coorPoint.reserve( 3 );
    coorPoint.push_back( 0.0 );
    coorPoint.push_back( 0.0 );
    coorPoint.push_back( 0.0 );

    ASTERINTEGER nbPoinQuad = 0;
    VectorReal poinQuadReal;
    poinQuadReal.reserve( 3 * 56 );
    CALLO_QUADPOINCOORWRAP( _mesh->getName(), _currentCoordinates->getName(), getBasename(),
                            &indexPair, &nbPoinQuad, poinQuadReal.data() );

    if ( nbPoinQuad != 0 ) {
        returnValue.reserve( nbPoinQuad );
        for ( long int iPoinQuad = 0; iPoinQuad < nbPoinQuad; iPoinQuad++ ) {
            coorPoint[0] = poinQuadReal[3 * iPoinQuad];
            coorPoint[1] = poinQuadReal[3 * iPoinQuad + 1];
            coorPoint[2] = poinQuadReal[3 * iPoinQuad + 2];
            returnValue.push_back( coorPoint );
        }
    }

    return returnValue;
}

/** @brief Check for common nodes between slave and master side */
ASTERBOOL MeshPairing::hasCommonNodes() const {

    ASTERBOOL returnValue;
    auto mesh = getMesh();

    VectorLong commonNodes;
    VectorLong masterNodes = _masterNodes;
    VectorLong slaveNodes = _slaveNodes;

    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        VectorLong slaveNodes_gl;
        AsterMPI::all_gather( _slaveNodes, slaveNodes_gl );
        commonNodes = set_intersection( slaveNodes_gl, masterNodes );
#endif
    } else {
        commonNodes = set_intersection( slaveNodes, masterNodes );
    }

    ASTERINTEGER size_inter_gl = commonNodes.size();

    // share error
#ifdef ASTER_HAVE_MPI
    if ( mesh->isParallel() ) {
        ASTERINTEGER size_inter_lc = size_inter_gl;
        size_inter_gl = AsterMPI::sum( size_inter_lc );
    }
#endif
    returnValue = false;
    if ( size_inter_gl > 0 ) {
        returnValue = true;
    }

    return returnValue;
}

VectorLong MeshPairing::getMasterCellsFromNode( const ASTERINTEGER &nodeIndex ) const {
    auto vct = ( *_masterInverseConnectivity )[nodeIndex + 1]->toVector();
    std::transform( vct.begin(), vct.end(), vct.begin(), [this]( ASTERINTEGER k ) -> ASTERINTEGER {
        return k > 0 ? _masterCells[k - 1] : 0;
    } );
    return vct;
}

VectorLong MeshPairing::getSlaveCellsFromNode( const ASTERINTEGER &nodeIndex ) const {
    auto vct = ( *_slaveInverseConnectivity )[nodeIndex + 1]->toVector();
    std::transform( vct.begin(), vct.end(), vct.begin(), [this]( ASTERINTEGER k ) -> ASTERINTEGER {
        return k > 0 ? _slaveCells[k - 1] : 0;
    } );
    return vct;
}

VectorLong MeshPairing::getMasterCellNeighbors( const ASTERINTEGER &cellIndex ) const {
    ASTERINTEGER ind_min = *std::min_element( _masterCells.begin(), _masterCells.end() );
    ASTERINTEGER ind_max = *std::max_element( _masterCells.begin(), _masterCells.end() );

    if ( cellIndex < ind_min || cellIndex > ind_max )
        throw std::out_of_range( " the master cell's number should be"
                                 " between " +
                                 std::to_string( ind_min ) + " and " + std::to_string( ind_max ) );

    auto vct = ( *_masterNeighbors )[cellIndex - ind_min + 1]->toVector();
    vct.erase( std::remove_if( vct.begin(), vct.end(),
                               []( ASTERINTEGER &cellIndex ) { return cellIndex == 0; } ),
               vct.end() );
    std::transform( vct.begin(), vct.end(), vct.begin(),
                    []( ASTERINTEGER k ) -> ASTERINTEGER { return k - 1; } );
    return vct;
}

VectorLong MeshPairing::getSlaveCellNeighbors( const ASTERINTEGER &cellIndex ) const {

    ASTERINTEGER ind_min = *std::min_element( _slaveCells.begin(), _slaveCells.end() );
    ASTERINTEGER ind_max = *std::max_element( _slaveCells.begin(), _slaveCells.end() );

    if ( cellIndex < ind_min || cellIndex > ind_max )
        throw std::out_of_range( " the slave cell's number should be"
                                 " between " +
                                 std::to_string( ind_min ) + " and " + std::to_string( ind_max ) );

    auto vct = ( *_slaveNeighbors )[cellIndex - ind_min + 1]->toVector();
    vct.erase( std::remove_if( vct.begin(), vct.end(),
                               []( ASTERINTEGER &cellIndex ) { return cellIndex == 0; } ),
               vct.end() );
    std::transform( vct.begin(), vct.end(), vct.begin(),
                    []( ASTERINTEGER k ) -> ASTERINTEGER { return k - 1; } );
    return vct;
}

ASTERBOOL MeshPairing::buildSlaveCellsVolu() {
    auto mesh = getMesh();
    auto nbCells = mesh->getNumberOfCells();

    auto invCon = mesh->getInverseConnectivity();
    invCon->build();

    auto exp = mesh->getConnectivityExplorer();

    for ( auto &cellId : _slaveCells ) {
        const auto cell = exp[cellId];
        VectorLong candidat;
        for ( const auto nodeId : cell ) {
            auto listCells = ( *invCon )[nodeId + 1]->toVector();

            std::for_each( listCells.begin(), listCells.end(), []( ASTERINTEGER &d ) { d -= 1; } );

            if ( candidat.empty() ) {
                candidat = listCells;
            } else {
                auto tmp = set_intersection( candidat, listCells );
                AS_ASSERT( !tmp.empty() );
                candidat = tmp;
            }
        }

        if ( candidat.size() == 2 ) {
            ASTERINTEGER cellVolu;
            if ( candidat[0] == cellId ) {
                cellVolu = candidat[1];
            } else {
                cellVolu = candidat[0];
            }

            _slavSurf2Volu[cellId] = cellVolu;
        }
    }

    return true;
};

void MeshPairing::setExcludedSlaveGroupOfCells( const VectorString &groupsName ) {
    _excludedSlaveCells = groupsName;
};

void MeshPairing::setExcludedSlaveGroupOfNodes( const VectorString &groupsName ) {
    _excludedSlaveNodes = groupsName;
};

ASTERBOOL MeshPairing::build() {
    if ( getMesh()->hasGroupOfCells( _slaveCellsGroup ) ) {
        _slaveNodes = getMesh()->getNodesFromCells( _slaveCellsGroup, false, true );
        _slaveCells = getMesh()->getCells( _slaveCellsGroup );
    } else {
        throw std::runtime_error( "The given group " + _slaveCellsGroup +
                                  " doesn't exist in mesh" );
    }
    if ( getMesh()->hasGroupOfCells( _masterCellsGroup ) ) {
        _masterNodes = getMesh()->getNodesFromCells( _masterCellsGroup, false, true );
        _masterCells = getMesh()->getCells( _masterCellsGroup );
    } else {
        throw std::runtime_error( "The given group " + _masterCellsGroup +
                                  " doesn't exist in mesh" );
    }

    if ( _excludedSlaveCells.size() != 0 ) {
        for ( auto &groupName : _excludedSlaveCells ) {
            if ( !( getMesh()->hasGroupOfCells( groupName ) ) ) {
                throw std::runtime_error( "The group " + groupName + " doesn't exist in mesh" );
            }

            VectorLong sans_gr_i = _mesh->getCells( groupName );
            auto it = _slaveCellsExcluded.end();
            _slaveCellsExcluded.insert( it, sans_gr_i.begin(), sans_gr_i.end() );
        }

        _slaveCells = set_difference( _slaveCells, _slaveCellsExcluded );
        AS_ASSERT( _slaveCells.size() > 0 );
        _slaveNodes = getMesh()->getNodesFromCells( _slaveCells );
    }

    if ( _excludedSlaveNodes.size() != 0 ) {
        // Build inverse connectivity
        AS_ASSERT( buildInverseConnectivity() );

        for ( auto &groupName : _excludedSlaveNodes ) {
            if ( !( getMesh()->hasGroupOfNodes( groupName ) ) ) {
                throw std::runtime_error( "The group " + groupName + " doesn't exist in mesh" );
            }

            VectorLong sans_gr_i = _mesh->getNodes( groupName );
            for ( auto &node : sans_gr_i ) {
                auto cellsToExclude = this->getSlaveCellsFromNode( node );
                auto it = _slaveCellsExcluded.end();
                _slaveCellsExcluded.insert( it, cellsToExclude.begin(), cellsToExclude.end() );
            }
        }

        _slaveCells = set_difference( _slaveCells, _slaveCellsExcluded );
        AS_ASSERT( _slaveCells.size() > 0 );
        _slaveNodes = getMesh()->getNodesFromCells( _slaveCells );
    }

    _zoneHaveBeenDefined = true;

    // Build inverse connectivity
    AS_ASSERT( buildInverseConnectivity() );

    // Build master and slave cells neighbors
    AS_ASSERT( buildCellsNeighbors() );

    // Build surface to volume slave cell mapping
    AS_ASSERT( buildSlaveCellsVolu() );

    return true;
}

void MeshPairing::checkNormals( const ModelPtr _model ) const {
    if ( _zoneHaveBeenDefined ) {
        if ( _model->getMesh() == getMesh() ) {
            CALL_CHECKNORMALS( _model->getName().c_str(),
                               ljust( getSlaveGroupOfCells(), 24, ' ' ).c_str(),
                               ljust( getMasterGroupOfCells(), 24, ' ' ).c_str() );
        } else {
            throw std::runtime_error( "Mesh is different in model" );
        }

    } else {
        throw std::runtime_error( "Zone is undefined" );
    }
}
