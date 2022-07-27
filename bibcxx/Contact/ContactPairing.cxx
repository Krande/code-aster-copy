/**
 * @file ContactPairing.cxx
 * @brief Implementation de Contact
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

ContactPairing::ContactPairing( const std::string name, const ContactNewPtr cont )
    : DataStructure( name, 8, "PAIRING_SD" ), _contDefi( cont ), _mesh( cont->getMesh() ) {
    if ( !_mesh || _mesh->isParallel() )
        raiseAsterError( "Mesh is empty or is parallel " );

    _newCoordinates = std::make_shared< MeshCoordinatesField >( *( _mesh->getCoordinates() ) );

    // be sure that zones is not empty and get size of zones and resize
    int size_zones = _contDefi->getNumberOfContactZones();
    if ( size_zones == 0 )
        raiseAsterError( "ContactZone vector is empty " );

    // resize pairing quantities
    _nbPairs.resize( size_zones );
    _listOfPairs.resize( size_zones );
    _nbIntersectionPoints.resize( size_zones );
    _slaveIntersectionPoints.resize( size_zones );
};

ASTERBOOL ContactPairing::computeZone( ASTERINTEGER i ) {

    if ( i < 0 || i >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }

    CALL_JEMARQ();

    auto zone = _contDefi->getContactZone( i );
    AS_ASSERT( !zone->hasSmoothing() );

    // get and define some input parameters
    VectorLong eleMaster = zone->getMasterCells();
    VectorLong NodesMaster = zone->getMasterNodes();
    VectorLong eleSlave = zone->getSlaveCells();
    ASTERINTEGER nbCellMaster = eleMaster.size();
    ASTERINTEGER nbNodeMaster = NodesMaster.size();
    ASTERINTEGER nbCellSlave = eleSlave.size();
    std::string pair_method;

    // update the numbering for fortran
    std::for_each( eleMaster.begin(), eleMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( NodesMaster.begin(), NodesMaster.end(), []( ASTERINTEGER &d ) { d += 1; } );
    std::for_each( eleSlave.begin(), eleSlave.end(), []( ASTERINTEGER &d ) { d += 1; } );

    // get pairing method
    auto variant = zone->getPairingParameter()->getAlgorithm();
    if ( variant == PairingAlgo::Mortar ) {
        pair_method = ljust( "RAPIDE", 24, ' ' );
    } else {
        AS_ABORT( "Not expected" );
    }

    // tolerence
    ASTERDOUBLE pair_tole = 1e-8;

    // set pairs numbers to 0
    ASTERINTEGER nb_pairs = 0;

    // output paramaters as C pointers
    auto pairs = JeveuxVectorLong( "&&LISTPAIRS" );
    auto nbInterPoints = JeveuxVectorLong( "&&NBPAIRS" );
    auto interSlavePoints = JeveuxVectorReal( "&&INTERSLPTS" );

    CALLO_APLCPGN( _mesh->getName(), _newCoordinates->getName(), zone->getName(), pair_method,
                   &pair_tole, &nbCellMaster, eleMaster.data(), &nbCellSlave, eleSlave.data(),
                   NodesMaster.data(), &nbNodeMaster, &nb_pairs, ljust( pairs->getName(), 19, ' ' ),
                   ljust( nbInterPoints->getName(), 19, ' ' ),
                   ljust( interSlavePoints->getName(), 19, ' ' ) );

    // clearZone
    this->clearZone( i );

    // fill the pairing quantities
    _nbPairs[i] = nb_pairs;
    _listOfPairs[i] = pairs->toVector();
    _nbIntersectionPoints[i] = nbInterPoints->toVector();
    _slaveIntersectionPoints[i] = interSlavePoints->toVector();

    // update numerotation

    std::transform( _listOfPairs[i].begin(), _listOfPairs[i].end(), _listOfPairs[i].begin(),
                    []( ASTERINTEGER &i ) -> ASTERINTEGER { return --i; } );

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

void ContactPairing::clearZone( ASTERINTEGER i ) {

    // swap is recommended to release memory
    _listOfPairs[i].clear();
    _nbIntersectionPoints[i].clear();
    _slaveIntersectionPoints[i].clear();

    _nbPairs.at( i ) = 0;
}

std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >
ContactPairing::getListOfPairsOfZone( ASTERINTEGER zone_index ) const {

    std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > > tmp;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( zone_index );
    tmp.reserve( nbPairs );

    for ( auto i = 0; i < nbPairs; i++ ) {
        tmp.push_back( std::make_pair( _listOfPairs[zone_index][2 * i],
                                       _listOfPairs[zone_index][2 * i + 1] ) );
    }
    return tmp;
}

std::vector< VectorReal >
ContactPairing::getSlaveIntersectionPoints( ASTERINTEGER zone_index ) const {

    std::vector< VectorReal > ret;
    ASTERINTEGER nbPairs = getNumberOfPairsOfZone( zone_index );
    ret.reserve( nbPairs );

    auto iter = _slaveIntersectionPoints[zone_index].begin();
    for ( auto i = 0; i < nbPairs; i++ ) {
        ret.push_back( VectorReal( iter + 16 * i, iter + 16 * ( i + 1 ) ) );
    }

    return ret;
}

void ContactPairing::buildFiniteElementDescriptor() {

    CALL_JEMARQ();

    auto model = _contDefi->getModel();
    auto mesh = getMesh();

    _fed = std::make_shared< FiniteElementDescriptor >( mesh );
    _fed->setModel( model );

    const ASTERINTEGER nbZoneCont = _contDefi->getNumberOfContactZones();
    const ASTERINTEGER nbContPairTot = this->getNumberOfPairs();

    using VectorPairLong = std::vector< std::pair< ASTERINTEGER, ASTERINTEGER > >;

    ASTERINTEGER nbType = 33;

    std::vector< std::array< ASTERINTEGER, 5 > > listType;
    listType.assign( nbType, std::array< ASTERINTEGER, 5 >( { 0, 0, 0, 0, 0 } ) );

    std::vector< VectorLong > listNodes;
    listNodes.reserve( nbContPairTot );

    std::vector< VectorPairLong > listContElem;
    listContElem.reserve( 2 * nbZoneCont );

    ASTERINTEGER modelDim = model->getGeometricDimension();
    ASTERLOGICAL lAxis = model->existsAxis();

    auto meshConnectivty = mesh->getConnectivity();

    SetLong slaveNodePaired, slaveCellPaired;
    ASTERINTEGER iContPair = 0;

    // Create virtual cells for pairing cells
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        auto iZonePairing = this->getListOfPairsOfZone( iZone );
        auto nbContPairZone = this->getNumberOfPairsOfZone( iZone );

        auto zone = _contDefi->getContactZone( iZone );
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFrot = zone->getFrictionParameter()->hasFriction();

        VectorPairLong listContElemZone;
        listContElemZone.reserve( nbContPairZone );

        /*loop on  pair of iZone*/
        for ( int iPair = 0; iPair < nbContPairZone; iPair++ ) {

            /*get pairing of current zone*/
            auto [slavCellNume, mastCellNume] = iZonePairing[iPair];
            slaveCellPaired.insert( slavCellNume );

            /*get slave and master geom type*/
            auto typgSlavName = ljust( mesh->getCellTypeName( slavCellNume ), 8, ' ' );
            auto typgMastName = ljust( mesh->getCellTypeName( mastCellNume ), 8, ' ' );

            /*call mmelemdata_c*/
            ASTERINTEGER typgContNume = 0, typfContNume = 0, typfFrotNume = 0, typeElem;
            ASTERINTEGER nbNodesCell = 0, elemIndx = 0;

            if ( contAlgo == ContactAlgo::Lagrangian ) {
                CALLO_MMELEM_DATA_LAGA( &lAxis, typgSlavName, typgMastName, &nbType, &nbNodesCell,
                                        &typgContNume, &typfContNume, &typfFrotNume, &elemIndx );
            } else {
                AS_ABORT( "Not implemented" );
            }

            auto &info = listType.at( elemIndx - 1 );

            if ( lFrot ) {
                info[1] += 1;
                info[3] = typfFrotNume;
                typeElem = typfFrotNume;
            } else {
                info[0] += 1;
                typeElem = typfContNume;
            }

            info[2] = typfContNume;
            info[4] = typgContNume;

            listContElemZone.push_back( std::make_pair( typeElem, ++iContPair ) );

            /* get nodes (be carefull with +1 ) */
            auto slav_cell_con = ( *meshConnectivty )[slavCellNume + 1];
            auto mast_cell_con = ( *meshConnectivty )[mastCellNume + 1];

            VectorLong toCopy;
            toCopy.reserve( slav_cell_con.size() + mast_cell_con.size() + 1 );

            /*Copy slave nodes*/
            auto toAdd = slav_cell_con.toVector();
            toCopy.insert( toCopy.end(), toAdd.begin(), toAdd.end() );
            slaveNodePaired.insert( toAdd.begin(), toAdd.end() );

            /*Copy master nodes*/
            toAdd = mast_cell_con.toVector();
            toCopy.insert( toCopy.end(), toAdd.begin(), toAdd.end() );

            /*Copy contact element type*/
            toCopy.push_back( typeElem );

            listNodes.push_back( toCopy );
        }
        listContElem.push_back( listContElemZone );
    }

    // Create Virtual cell for orphelan nodes
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        auto iZonePairing = this->getListOfPairsOfZone( iZone );
        auto nbContPairZone = this->getNumberOfPairsOfZone( iZone );

        auto zone = _contDefi->getContactZone( iZone );
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFrot = zone->getFrictionParameter()->hasFriction();

        auto slaveCells = zone->getSlaveCells();

        VectorPairLong listContElemZone;
        listContElemZone.reserve( slaveCells.size() );

        // Find slave cells that are not paired
        for ( auto &slavCellNume : slaveCells ) {
            if ( slaveCellPaired.count( slavCellNume ) == 0 ) {
                slaveCellPaired.insert( slavCellNume );
                auto slav_cell_con = ( *meshConnectivty )[slavCellNume + 1].toVector();
                auto cellType = trim( mesh->getCellTypeName( slavCellNume ) );
                ASTERINTEGER nno_lgar = 0, nno = 0;

                if ( cellType == "SEG2" ) {
                    nno_lgar = 2;
                } else if ( cellType == "SEG3" ) {
                    nno_lgar = 3;
                } else if ( cellType == "TRIA3" || cellType == "TRIA6" || cellType == "TRIA7" ) {
                    nno_lgar = 3;
                } else if ( cellType == "QUAD4" || cellType == "QUAD8" || cellType == "QUAD9" ) {
                    nno_lgar = 4;
                } else {
                    AS_ABORT( cellType + " not supported" );
                }

                for ( auto &nodeNume : slav_cell_con ) {
                    nno++;
                    if ( slaveNodePaired.count( nodeNume ) == 0 ) {
                        slaveNodePaired.insert( nodeNume );

                        auto typgSlavName = ljust( "POI1", 8, ' ' );
                        std::string typgMastName;
                        if ( nno <= nno_lgar ) {
                            typgMastName = ljust( "LAG" + std::to_string( modelDim ), 8, ' ' );
                        } else {
                            typgMastName = ljust( "NOLAG" + std::to_string( modelDim ), 8, ' ' );
                        }

                        /*call mmelemdata_c*/
                        ASTERINTEGER typgContNume = 0, typfContNume = 0, typfFrotNume = 0;
                        ASTERINTEGER nbNodesCell = 0, elemIndx = 0, typeElem;

                        if ( contAlgo == ContactAlgo::Lagrangian ) {
                            CALLO_MMELEM_DATA_LAGA( &lAxis, typgSlavName, typgMastName, &nbType,
                                                    &nbNodesCell, &typgContNume, &typfContNume,
                                                    &typfFrotNume, &elemIndx );
                        } else {
                            AS_ABORT( "Not implemented" );
                        }

                        AS_ASSERT( nbNodesCell == 1 );

                        auto &info = listType.at( elemIndx - 1 );

                        if ( lFrot ) {
                            info[1] += 1;
                            info[3] = typfFrotNume;
                            typeElem = typfFrotNume;
                        } else {
                            info[0] += 1;
                            typeElem = typfContNume;
                        }

                        info[2] = typfContNume;
                        info[4] = typgContNume;

                        listContElemZone.push_back( std::make_pair( typeElem, ++iContPair ) );
                        listNodes.push_back( VectorLong( { nodeNume, typeElem } ) );
                    }
                }
            }
        }
        if ( !listContElemZone.empty() ) {
            std::cout << "Not paired nodes: " << listContElemZone.size() << std::endl;
            listContElem.push_back( listContElemZone );
        }
    }
    slaveNodePaired.clear();
    slaveCellPaired.clear();

    /*FED building*/
    auto ContactResFEDNbno = _fed->getNumberOfVirtualNodesobj();
    ContactResFEDNbno->allocate( 1 );
    ( *ContactResFEDNbno )[0] = 0;

    /*NEMA building*/
    auto ContactResFEDNema = _fed->getNema();
    ContactResFEDNema->allocateContiguousNumbered( listNodes );
    listNodes.clear();

    /*LIEL building
    Size of LIEL object*/
    ASTERINTEGER nbGrel = 0, ligrcf_liel_lont = 0;
    for ( auto &info : listType ) {
        auto nbCont = info[0];
        auto nbFrot = info[1];

        ligrcf_liel_lont += nbCont + nbFrot;

        if ( nbCont > 0 ) {
            nbGrel += 1;
        }
        if ( nbFrot > 0 ) {
            nbGrel += 1;
        }
    }
    ligrcf_liel_lont += nbGrel;

    /*Create LIEL object*/
    auto ContactResFEDLiel = _fed->getListOfGroupOfElements();
    ContactResFEDLiel->allocateContiguousNumbered( nbGrel, ligrcf_liel_lont, Variable );

    _pair2Zone.clear();

    for ( auto &info : listType ) {
        auto nbCont = info[0];
        auto nbFrot = info[1];

        /* Create and add Contact element*/
        if ( nbCont > 0 ) {
            auto typfContNume = info[2];

            VectorLong toCopy;
            toCopy.reserve( nbCont );

            ASTERINTEGER iZone = 0;
            for ( auto &listContElemZone : listContElem ) {
                /*loop on  pair of iZone*/
                for ( auto &[typContNume, iContPair] : listContElemZone ) {
                    if ( typContNume == typfContNume ) {
                        toCopy.push_back( -iContPair );
                        _pair2Zone[iContPair - 1] = iZone;
                    }
                }
                iZone++;
            }
            toCopy.push_back( typfContNume );
            ContactResFEDLiel->push_back( toCopy );
        }

        /* Create friction element*/
        if ( nbFrot > 0 ) {
            auto typfFrotNume = info[3];

            VectorLong toCopy;
            toCopy.reserve( nbCont );

            ASTERINTEGER iZone = 0;
            for ( auto &listContElemZone : listContElem ) {
                /*loop on  pair of iZone*/
                for ( auto &[typFrotNume, iContPair] : listContElemZone ) {
                    if ( typFrotNume == typfFrotNume ) {
                        toCopy.push_back( -iContPair );
                        _pair2Zone[iContPair - 1] = iZone;
                    }
                }
                iZone++;
            }
            toCopy.push_back( typfFrotNume );
            ContactResFEDLiel->push_back( toCopy );
        }
    }

    /*Create LGRF object*/
    auto paramToCopy = model->getFiniteElementDescriptor()->getParameters();
    paramToCopy->updateValuePointer();
    auto parameters = _fed->getParameters();
    parameters->allocate( 3 );

    ( *parameters )[0] = mesh->getName();
    ( *parameters )[1] = model->getName();
    ( *parameters )[2] = ( *paramToCopy )[2];
    auto docu = paramToCopy->getInformationParameter();
    parameters->setInformationParameter( docu );

    CALLO_ADALIG_WRAP( _fed->getName() );
    bool l_calc_rigi = false;
    CALLO_INITEL( _fed->getName(), (ASTERLOGICAL *)&l_calc_rigi );

    CALL_JEDEMA();
};
