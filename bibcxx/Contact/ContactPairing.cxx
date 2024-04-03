/**
 * @file ContactPairing.cxx
 * @brief Implementation de Contact
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "Contact/ContactPairing.h"

#include "Messages/Messages.h"

void ContactPairing::resizePairing( const int nbZoneCont ) {
    if ( nbZoneCont == 0 )
        raiseAsterError( "ContactZone vector is empty " );
    _nbPairs.resize( nbZoneCont );
    _listOfPairs.resize( nbZoneCont );
    _nbIntersectionPoints.resize( nbZoneCont );
    _slaveIntersectionPoints.resize( nbZoneCont );
}

ContactPairing::ContactPairing( const std::string name, const ContactNewPtr cont )
    : DataStructure( name, 8, "PAIRING_SD" ), _contDefi( cont ), _mesh( cont->getMesh() ) {
    if ( !_mesh || _mesh->isParallel() )
        raiseAsterError( "Mesh is empty or is parallel " );

    _currentCoordinates = std::make_shared< MeshCoordinatesField >( *( _mesh->getCoordinates() ) );

    // be sure that zones is not empty and get size of zones and resize
    int nbZoneCont = _contDefi->getNumberOfContactZones();
    if ( nbZoneCont == 0 )
        raiseAsterError( "ContactZone vector is empty " );

    // Resize pairing quantities
    resizePairing( nbZoneCont );
};

ASTERBOOL ContactPairing::computeZone( ASTERINTEGER indexZone ) {

    if ( indexZone < 0 || indexZone >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }

    CALL_JEMARQ();

    auto zone = _contDefi->getContactZone( indexZone );
    AS_ASSERT( !zone->hasSmoothing() );

    // Get and define some input parameters
    VectorLong masterCells = zone->getMasterCells();
    VectorLong masterNodes = zone->getMasterNodes();
    VectorLong slaveCells = zone->getSlaveCells();
    ASTERINTEGER nbCellMaster = masterCells.size();
    ASTERINTEGER nbNodeMaster = masterNodes.size();
    ASTERINTEGER nbCellSlave = slaveCells.size();
    std::string pair_method;

    // Update the numbering for fortran
    std::for_each( masterCells.begin(), masterCells.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( masterNodes.begin(), masterNodes.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( slaveCells.begin(), slaveCells.end(), []( ASTERINTEGER &d ) { d += 1; } );

    // Get pairing method
    auto variant = zone->getPairingParameter()->getAlgorithm();
    if ( variant == PairingAlgo::Mortar ) {
        pair_method = ljust( "RAPIDE", 24, ' ' );
    } else {
        AS_ABORT( "Not expected" );
    }

    // Get distance ratio
    auto dist_pairing = zone->getPairingParameter()->getDistanceRatio();

    // Tolerance for pairing
    ASTERDOUBLE pair_tole = 1e-8;

    // Set pairs numbers to 0
    ASTERINTEGER nb_pairs = 0;

    // Main routine for pairing
    auto pairs = JeveuxVectorLong( "&&LISTPAIRS" );
    auto nbInterPoints = JeveuxVectorLong( "&&NBPAIRS" );
    auto interSlavePoints = JeveuxVectorReal( "&&INTERSLPTS" );

    CALLO_APLCPGN( _mesh->getName(), _currentCoordinates->getName(), zone->getName(), pair_method,
                   &pair_tole, &dist_pairing, &nbCellMaster, masterCells.data(), &nbCellSlave,
                   slaveCells.data(), masterNodes.data(), &nbNodeMaster, &nb_pairs,
                   ljust( pairs->getName(), 19, ' ' ), ljust( nbInterPoints->getName(), 19, ' ' ),
                   ljust( interSlavePoints->getName(), 19, ' ' ) );

    // clearZone
    this->clearZone( indexZone );

    // fill the pairing quantities
    _nbPairs[indexZone] = nb_pairs;
    _listOfPairs[indexZone] = pairs->toVector();
    _nbIntersectionPoints[indexZone] = nbInterPoints->toVector();
    _slaveIntersectionPoints[indexZone] = interSlavePoints->toVector();

    // update numerotation
    std::transform( _listOfPairs[indexZone].begin(), _listOfPairs[indexZone].end(),
                    _listOfPairs[indexZone].begin(),
                    []( ASTERINTEGER &indexZone ) -> ASTERINTEGER { return --indexZone; } );

    CALL_JEDEMA();

    return true;
}

ASTERBOOL ContactPairing::compute() {
    // Pairing
    for ( int i = 0; i < _contDefi->getNumberOfContactZones(); i++ ) {
        computeZone( i );
    }

    // Build FED
    this->buildFiniteElementDescriptor();

    return true;
}

void ContactPairing::clearZone( ASTERINTEGER indexZone ) {

    _listOfPairs[indexZone].clear();
    _nbIntersectionPoints[indexZone].clear();
    _slaveIntersectionPoints[indexZone].clear();
    _nbPairs.at( indexZone ) = 0;
}

std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >
ContactPairing::getListOfPairsOfZone( ASTERINTEGER indexZone ) const {

    std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > tmp;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( indexZone );
    tmp.reserve( nbPairs );

    for ( auto iPair = 0; iPair < nbPairs; iPair++ ) {
        tmp.push_back( std::make_pair( _listOfPairs[indexZone][2 * iPair],
                                       _listOfPairs[indexZone][2 * iPair + 1] ) );
    }

    return tmp;
}

std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > ContactPairing::getListOfPairs() const {

    std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > tmp;
    ASTERINTEGER nbPairs = getNumberOfPairs();
    tmp.reserve( nbPairs );

    for ( int indexZone = 0; indexZone < _contDefi->getNumberOfContactZones(); indexZone++ ) {
        auto nbPairs = getNumberOfPairsOfZone( indexZone );

        for ( auto iPair = 0; iPair < nbPairs; iPair++ ) {
            tmp.push_back( std::make_pair( _listOfPairs[indexZone][2 * iPair],
                                           _listOfPairs[indexZone][2 * iPair + 1] ) );
        }
    }

    return tmp;
}

std::vector< VectorReal >
ContactPairing::getSlaveIntersectionPoints( ASTERINTEGER indexZone ) const {

    std::vector< VectorReal > ret;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( indexZone );
    ret.reserve( nbPairs );

    auto iter = _slaveIntersectionPoints[indexZone].begin();
    for ( auto i = 0; i < nbPairs; i++ ) {
        ret.push_back( VectorReal( iter + 16 * i, iter + 16 * ( i + 1 ) ) );
    }

    return ret;
}

ASTERINTEGER ContactPairing::getContCellIndx( const ContactAlgo contAlgo,
                                              std::string slavCellTypeName,
                                              std::string mastCellTypeName ) {
    ASTERINTEGER cellIndx = -1;

    if ( contAlgo == ContactAlgo::Lagrangian ) {
        for ( int iContType = 0; iContType < contLagrType; iContType++ ) {
            if ( slavCellTypeName == contCellLagr[iContType].slavCellType ) {
                if ( mastCellTypeName == contCellLagr[iContType].mastCellType ) {
                    AS_ASSERT( cellIndx == -1 )
                    cellIndx = iContType;
                }
            }
        }
    } else if ( contAlgo == ContactAlgo::Nitsche ) {
        for ( int iContType = 0; iContType < contNitsType; iContType++ ) {
            if ( slavCellTypeName == contCellNits[iContType].slavCellType ) {
                if ( mastCellTypeName == contCellNits[iContType].mastCellType ) {
                    AS_ASSERT( cellIndx == -1 )
                    cellIndx = iContType;
                }
            }
        }
    } else {
        AS_ABORT( "Not implemented" );
    };

    return cellIndx;
}

ASTERINTEGER ContactPairing::getContCellType( const ContactAlgo contAlgo,
                                              const ASTERINTEGER cellIndx, const bool lAxis,
                                              const bool lFric ) {

    AS_ASSERT( cellIndx != -1 );

    ASTERINTEGER contTypeNume = -1;

    // Get finite element descriptor of model
    auto modelFEDesc = _contDefi->getModel()->getFiniteElementDescriptor();

    // Get name of type of contact cell
    std::string contTypeName;
    if ( contAlgo == ContactAlgo::Lagrangian ) {
        if ( lFric ) {
            contTypeName = contCellLagr[cellIndx].fricElemType;
        } else {
            contTypeName = contCellLagr[cellIndx].contElemType;
        }
    } else if ( contAlgo == ContactAlgo::Nitsche ) {
        if ( lFric ) {
            contTypeName = contCellNits[cellIndx].fricElemType;
        } else {
            contTypeName = contCellNits[cellIndx].contElemType;
        }
    } else {
        AS_ABORT( "Not implemented" );
    };

    if ( lAxis ) {
        contTypeName.append( "A" );
    }

    // Get index of type of contact cell
    contTypeNume = modelFEDesc->getElemTypeNume( contTypeName );
    return contTypeNume;
}

void ContactPairing::createVirtualElemForContact( const ASTERLOGICAL lAxis, const int nbZoneCont,
                                                  MapLong &contactElemType,
                                                  const JeveuxCollectionLong meshConnectivity,
                                                  std::vector< VectorLong > &listContElem,
                                                  std::vector< VectorPairLong > &listContType,
                                                  ASTERINTEGER &iContPair, SetLong &slaveNodePaired,
                                                  SetLong &slaveCellPaired ) {

    auto mesh = getMesh();

    // Loop on contact zones
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        // Get current zone
        auto zone = _contDefi->getContactZone( iZone );

        // Get contact parameters for this zone
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFric = zone->getFrictionParameter()->hasFriction();

        // Get pairing for this zone
        auto surf2Volu = zone->getSlaveCellsSurfToVolu();
        auto iZonePairing = this->getListOfPairsOfZone( iZone );
        auto nbContPairZone = this->getNumberOfPairsOfZone( iZone );

        // Create vector of (virtual) contact cells for this zone
        VectorPairLong listContTypeZone;
        listContTypeZone.reserve( nbContPairZone );

        // Loop on pairs in zone
        for ( int iPair = 0; iPair < nbContPairZone; iPair++ ) {

            // Get pairing of current zone
            auto [slavCellNume, mastCellNume] = iZonePairing[iPair];

            // Get cell slave to construct (is volumic cell for Nitsche)
            auto slavCellUsedNume = slavCellNume;
            if ( contAlgo == ContactAlgo::Nitsche ) {
                slavCellUsedNume = surf2Volu[slavCellNume];
            }

            // Get slave and master cell type
            auto slavCellTypeName = mesh->getCellTypeName( slavCellUsedNume );
            auto mastCellTypeName = mesh->getCellTypeName( mastCellNume );

            // Get index of contact cell (in fixed lists)
            ASTERINTEGER cellIndx = getContCellIndx( contAlgo, slavCellTypeName, mastCellTypeName );
            AS_ASSERT( cellIndx != -1 );

            // Get index of type of contact cell
            ASTERINTEGER typeElemNume = getContCellType( contAlgo, cellIndx, lAxis, lFric );
            AS_ASSERT( typeElemNume != -1 );

            // Add contact element to list
            if ( contactElemType.count( typeElemNume ) == 0 ) {
                contactElemType[typeElemNume] = 0;
            }
            contactElemType[typeElemNume] += 1;

            listContTypeZone.push_back( std::make_pair( typeElemNume, ++iContPair ) );

            // Get nodes
            auto slav_cell_con = ( *meshConnectivity )[slavCellUsedNume + 1];
            auto mast_cell_con = ( *meshConnectivity )[mastCellNume + 1];

            // Contact element on zone
            VectorLong contactElemZone;
            contactElemZone.reserve( slav_cell_con->size() + mast_cell_con->size() + 1 );

            // Copy slave nodes to contact element
            auto toAdd = slav_cell_con->toVector();
            contactElemZone.insert( contactElemZone.end(), toAdd.begin(), toAdd.end() );

            // Add slave nodes to list of paired nodes
            slaveNodePaired.insert( toAdd.begin(), toAdd.end() );

            // Add slave cell to list of paired cells
            slaveCellPaired.insert( slavCellNume );

            // Copy master nodes to contact element
            toAdd = mast_cell_con->toVector();
            contactElemZone.insert( contactElemZone.end(), toAdd.begin(), toAdd.end() );

            // Add type of contact element
            contactElemZone.push_back( typeElemNume );

            // Add contact element to all contact elements
            listContElem.push_back( contactElemZone );
        }
        // Add contact elements of zone
        listContType.push_back( listContTypeZone );
    }
}

void ContactPairing::createVirtualElemForOrphelanNodes(
    const ASTERLOGICAL lAxis, const int nbZoneCont, MapLong &contactElemType,
    const JeveuxCollectionLong meshConnectivity, std::vector< VectorLong > &listContElem,
    std::vector< VectorPairLong > &listContType, ASTERINTEGER &iContPair, SetLong &slaveNodePaired,
    SetLong &slaveCellPaired ) {

    // Get mesh
    auto mesh = getMesh();

    // Get model
    auto model = _contDefi->getModel();
    ASTERINTEGER modelDim = model->getGeometricDimension();

    // Loop on contact zones
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        // Get current zone
        auto zone = _contDefi->getContactZone( iZone );

        // Get contact parameters for this zone
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFric = zone->getFrictionParameter()->hasFriction();

        // Get pairing for this zone
        auto iZonePairing = this->getListOfPairsOfZone( iZone );
        auto nbContPairZone = this->getNumberOfPairsOfZone( iZone );

        // Get slave cells on this zone
        auto slaveCells = zone->getSlaveCells();

        // Create vector of (virtual) contact cells for this zone
        VectorPairLong listContTypeZone;
        listContTypeZone.reserve( slaveCells.size() );

        // Find slave cells that are not paired (because of LAGR_C)
        for ( auto &slavCellNume : slaveCells ) {
            if ( slaveCellPaired.count( slavCellNume ) == 0 ) {
                slaveCellPaired.insert( slavCellNume );

                // Get nodes of slave cell
                auto slav_cell_con = ( *meshConnectivity )[slavCellNume + 1]->toVector();

                // Get slave cell type
                auto slavCellTypeName = mesh->getCellTypeName( slavCellNume );

                // Select number of nodes with LAGR_C: P2/P1 for 3D and P2/P2 for 2D
                ASTERINTEGER nno_lgar = 0;

                if ( slavCellTypeName == "SEG2" ) {
                    nno_lgar = 2;
                } else if ( slavCellTypeName == "SEG3" ) {
                    nno_lgar = 3;
                } else if ( slavCellTypeName == "TRIA3" || slavCellTypeName == "TRIA6" ||
                            slavCellTypeName == "TRIA7" ) {
                    nno_lgar = 3;
                } else if ( slavCellTypeName == "QUAD4" || slavCellTypeName == "QUAD8" ||
                            slavCellTypeName == "QUAD9" ) {
                    nno_lgar = 4;
                } else {
                    AS_ABORT( slavCellTypeName + " not supported" );
                }

                // Loop on nodes on slave cell
                ASTERINTEGER nno = 0;
                for ( auto &nodeNume : slav_cell_con ) {
                    nno++;
                    // This node hasn't been paired
                    if ( slaveNodePaired.count( nodeNume ) == 0 ) {
                        slaveNodePaired.insert( nodeNume );

                        // Type of cell for slave: POI1
                        auto slavCellTypeName = "POI1";

                        // Type of cell for master
                        std::string mastCellTypeName;
                        if ( nno <= nno_lgar ) {
                            mastCellTypeName = "LAG" + std::to_string( modelDim );
                        } else {
                            mastCellTypeName = "NOLAG" + std::to_string( modelDim );
                        }

                        if ( contAlgo == ContactAlgo::Lagrangian ) {

                        } else if ( contAlgo == ContactAlgo::Nitsche ) {
                            continue;
                        } else {
                            AS_ABORT( "Not implemented" );
                        }

                        // Get index of contact cell (in fixed lists)
                        ASTERINTEGER cellIndx =
                            getContCellIndx( contAlgo, slavCellTypeName, mastCellTypeName );
                        AS_ASSERT( cellIndx != -1 );

                        // Number of nodes
                        ASTERINTEGER nbNodesCell = contCellLagr[cellIndx].nbNode;
                        AS_ASSERT( nbNodesCell == 1 );

                        // Get index of type of contact cell
                        ASTERINTEGER typeElemNume =
                            getContCellType( contAlgo, cellIndx, lAxis, lFric );
                        AS_ASSERT( typeElemNume != -1 )

                        // Add type of contact element
                        if ( contactElemType.count( typeElemNume ) == 0 ) {
                            contactElemType[typeElemNume] = 0;
                        }
                        contactElemType[typeElemNume] += 1;

                        // Add the new contact element
                        listContTypeZone.push_back( std::make_pair( typeElemNume, ++iContPair ) );

                        // Add the nodes of the new contact element
                        listContElem.push_back( VectorLong( { nodeNume, typeElemNume } ) );
                    }
                }
            }
        }
        if ( !listContTypeZone.empty() ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Not paired nodes: " << listContTypeZone.size() << std::endl;
#endif
            listContType.push_back( listContTypeZone );
        }
    }
};

void ContactPairing::buildFiniteElementDescriptor() {

    CALL_JEMARQ();

    // Model
    auto model = _contDefi->getModel();
    const ASTERLOGICAL lAxis = model->existsAxis();

    // Mesh
    auto mesh = getMesh();
    const JeveuxCollectionLong meshConnectivity = mesh->getConnectivity();

    // Get pairing parameters
    const ASTERINTEGER nbZoneCont = _contDefi->getNumberOfContactZones();
    const ASTERINTEGER nbContPairTot = this->getNumberOfPairs();

    // Number of cells for each type of contact cell
    MapLong contactElemType;

    // Create objets for nodes and cells
    VectorOfVectorsLong listContElem;
    listContElem.reserve( nbContPairTot );

    std::vector< VectorPairLong > listContType;
    listContType.reserve( 2 * nbZoneCont );

    // Objects
    SetLong slaveNodePaired, slaveCellPaired;

    // Index of current contact pair
    ASTERINTEGER iContPair = 0;

    // Create virtual elements for contact
    createVirtualElemForContact( lAxis, nbZoneCont, contactElemType, meshConnectivity, listContElem,
                                 listContType, iContPair, slaveNodePaired, slaveCellPaired );

    // Create virtual elements for orphelan nodes
    createVirtualElemForOrphelanNodes( lAxis, nbZoneCont, contactElemType, meshConnectivity,
                                       listContElem, listContType, iContPair, slaveNodePaired,
                                       slaveCellPaired );

    // Create finite element descriptor for virtual contact elements
    _fed = std::make_shared< FiniteElementDescriptor >( mesh );
    _fed->setModel( model );

    // Create list of virtual nodes (none !)
    _fed->setNumberOfVirtualNodes( 0 );

    // Create list of virtual elements (NEMA object)
    auto ContactResFEDNema = _fed->getVirtualCellsDescriptor();
    ContactResFEDNema->allocateContiguousNumbered( listContElem );

    // Number of groups of elements and length of FED for contact element
    ASTERINTEGER nbGrel = 0, ligrcf_liel_lont = 0;
    for ( auto &[type, size] : contactElemType ) {
        ligrcf_liel_lont += size;
        nbGrel += 1;
    }
    ligrcf_liel_lont += nbGrel;

    // Create list of elements (LIEL object)
    auto ContactResFEDLiel = _fed->getListOfGroupsOfElements();
    ContactResFEDLiel->allocateContiguousNumbered( nbGrel, ligrcf_liel_lont, Variable );

    // Clear map between zone and contact elements
    _pair2Zone.clear();

    // Add virtual element
    for ( auto &[type, size] : contactElemType ) {
        VectorLong contactElemZone;
        contactElemZone.reserve( size );

        ASTERINTEGER iZone = 0;
        for ( auto &listContTypeZone : listContType ) {
            for ( auto &[typeElemNume, iContPair] : listContTypeZone ) {
                if ( typeElemNume == type ) {
                    // Virtual cells = index of cells is negative
                    contactElemZone.push_back( -iContPair );
                    _pair2Zone[iContPair - 1] = iZone;
                }
            }
            iZone++;
        }
        contactElemZone.push_back( type );
        ContactResFEDLiel->push_back( contactElemZone );
    }

    // Get parameters from standard model
    auto paramToCopy = model->getFiniteElementDescriptor()->getParameters();
    paramToCopy->updateValuePointer();
    auto docu = paramToCopy->getInformationParameter();

    // Create LGRF object
    auto parameters = _fed->getParameters();
    parameters->allocate( 3 );
    ( *parameters )[0] = mesh->getName();
    ( *parameters )[1] = model->getName();
    ( *parameters )[2] = ( *paramToCopy )[2];
    parameters->setInformationParameter( docu );

    // Adapt FED
    CALLO_ADALIG_WRAP( _fed->getName() );
    bool l_calc_rigi = false;
    CALLO_INITEL( _fed->getName(), (ASTERLOGICAL *)&l_calc_rigi );

    // Final building
    _fed->build();

    // Clean
    listContElem.clear();
    slaveNodePaired.clear();
    slaveCellPaired.clear();

    CALL_JEDEMA();
};
