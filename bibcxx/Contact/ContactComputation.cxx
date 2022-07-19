
/**
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

#include "Contact/ContactComputation.h"

#include "DataFields/FieldOnNodes.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Messages/Messages.h"
#include "Numbering/DOFNumbering.h"

void ContactComputation::buildContactResFED( const ContactPairingPtr pairing ) {
    _fed = std::make_shared< FiniteElementDescriptor >( _contact->getMesh() );
    _fed->setModel( _contact->getModel() );

    const ASTERINTEGER nbZoneCont = _contact->getNumberOfContactZones();
    const ASTERINTEGER nbContPairTot = pairing->getNumberOfPairs();

    ASTERINTEGER nbType = 29;

    std::vector< std::array< ASTERINTEGER, 5 > > listType;
    listType.resize( nbType );

    std::vector< VectorLong > listNodes;
    listNodes.reserve( nbContPairTot );

    std::vector< VectorLong > listContElem;
    listContElem.reserve( nbZoneCont );

    auto mesh = _contact->getMesh();
    auto model = _contact->getModel();
    ASTERINTEGER modelDim = model->getGeometricDimension();
    ASTERLOGICAL lAxis = model->existsAxis();

    auto meshConnectivty = mesh->getConnectivity();

    /*loop on contact zone*/
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        auto iZonePairing = pairing->getListOfPairsOfZone( iZone );
        auto nbContPairZone = pairing->getNumberOfPairsOfZone( iZone );

        auto zone = _contact->getContactZone( iZone );
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFrot = zone->getFrictionParameter()->hasFriction();
        AS_ASSERT( !lFrot );

        VectorLong listContElemZone;
        listContElemZone.reserve( nbContPairZone );

        /*loop on  pair of iZone*/
        for ( int iPair = 0; iPair < nbContPairZone; iPair++ ) {

            /*get pairing of current zone*/
            auto [elemSlavNume, elemMastNume] = iZonePairing[iPair];

            /*get slave and master geom type*/
            auto typgSlavName = ljust( mesh->getCellTypeName( elemSlavNume ), 8, ' ' );
            auto typgMastName = ljust( mesh->getCellTypeName( elemMastNume ), 8, ' ' );

            /*call mmelemdata_c*/
            ASTERINTEGER typgContNume = 0, typfContNume = 0, typfFrotNume = 0;
            ASTERINTEGER nbNodeElem = 0, elemIndx = 0;

            if ( contAlgo == ContactAlgo::Lagrangian ) {
                CALLO_MMELEM_DATA_LAGA( &lAxis, &modelDim, typgSlavName, typgMastName, &nbType,
                                        &nbNodeElem, &typgContNume, &typfContNume, &typfFrotNume,
                                        &elemIndx );
            } else {
                AS_ABORT( "Not implemented" );
            }

            listContElemZone.push_back( typfContNume );

            auto &info = listType[elemIndx];

            if ( lFrot ) {
                info[1] += 1;
                info[3] = typfFrotNume;
            } else {
                info[0] += 1;
            }

            info[2] = typfContNume;
            info[4] = typgContNume;

            /* get nodes (be carefull with +1 ) */
            auto slav_cell_con = ( *meshConnectivty )[elemSlavNume + 1];
            auto mast_cell_con = ( *meshConnectivty )[elemMastNume + 1];

            VectorLong toCopy;
            toCopy.reserve( slav_cell_con.size() + mast_cell_con.size() + 1 );

            /*Copy slave nodes*/
            auto toAdd = slav_cell_con.toVector();
            toCopy.insert( toCopy.end(), toAdd.begin(), toAdd.end() );

            /*Copy master nodes*/
            toAdd = mast_cell_con.toVector();
            toCopy.insert( toCopy.end(), toAdd.begin(), toAdd.end() );

            /*Copy contact element type*/
            toCopy.push_back( typfContNume );

            listNodes.push_back( toCopy );
        }
        listContElem.push_back( listContElemZone );
    }

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
    auto ContactResFEDLiel = _fed->getListOfGroupOfCells();
    ContactResFEDLiel->allocateContiguousNumbered( nbGrel, ligrcf_liel_lont, Variable );

    for ( auto &info : listType ) {
        auto nbCont = info[0];
        auto nbFrot = info[1];

        /* Create and add Contact element*/
        if ( nbCont > 0 ) {
            auto typfContNume = info[2];

            VectorLong toCopy;
            ASTERINTEGER iContPair = 0;

            for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
                auto iZonePairing = pairing->getListOfPairsOfZone( iZone );
                auto &listContElemZone = listContElem[iZone];
                auto nbContPairZone = pairing->getNumberOfPairsOfZone( iZone );

                /*loop on  pair of iZone*/
                for ( auto &typContNume : listContElemZone ) {
                    iContPair += 1;
                    if ( typContNume == typfContNume ) {
                        toCopy.push_back( -iContPair );
                    }
                }
            }
            toCopy.push_back( typfContNume );
            ContactResFEDLiel->push_back( toCopy );
        }

        /* Create friction element*/
        if ( nbFrot != 0 ) {
            AS_ABORT( "Friction not implemented" );
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
}

std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr >
ContactComputation::geometricGap( const MeshCoordinatesFieldPtr coor ) const {

    // Prepare computing
    const std::string option( "GAP_GEOM" );
    CalculPtr _calcul = std::make_unique< Calcul >( option );

    _calcul->setFiniteElementDescriptor( _fed );

    // Add input fields
    _calcul->addInputField( "PGEOMER", coor );

    // Add output elementary
    auto gap_elem = std::make_shared< ElementaryVectorReal >();
    gap_elem->prepareCompute( option );

    auto igap_elem = std::make_shared< ElementaryVectorReal >();
    igap_elem->prepareCompute( option );

    _calcul->addOutputElementaryTerm( "PVECGAP", std::make_shared< ElementaryTermReal >() );
    _calcul->addOutputElementaryTerm( "PVEIGAP", std::make_shared< ElementaryTermReal >() );

    // compute
    _calcul->compute();

    // get elementary Term
    gap_elem->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVECGAP" ) );
    gap_elem->build();

    igap_elem->addElementaryTerm( _calcul->getOutputElementaryTerm( "PVEIGAP" ) );
    igap_elem->build();

    // Assemble
    auto gap = gap_elem->assemble();
    auto igap = igap_elem->assemble();

    // Divide by number of pairing
    AS_ASSERT( gap->isSimilarTo( *igap ) );
    auto nbDofs = gap->size();

    gap->updateValuePointers();
    igap->updateValuePointers();

    bool alarm = false;

    for ( ASTERINTEGER iDoF = 0; iDoF < nbDofs; iDoF++ ) {
        ASTERINTEGER indi = ASTERINTEGER( ( *igap )[iDoF] );

        if ( indi > 2 ) {
            alarm = true;
        }

        if ( indi == 0 ) {
            ( *gap )[iDoF] = std::nan( "" );
            ( *igap )[iDoF] = 0.0;
        } else {
            ( *gap )[iDoF] /= ASTERDOUBLE( indi );
            ( *igap )[iDoF] = 1.0;
        }
    }

    if ( alarm ) {
        UTMESS( "A", "CONTACT1_5" );
    }

    return std::make_pair( gap, igap );
}

/**
 * @brief Compute contact mortar matrix
 */
ElementaryMatrixDisplacementRealPtr ContactComputation::contactMortarMatrix() const {
    auto elemMatr = std::make_shared< ElementaryMatrixDisplacementReal >();

    const std::string option( "XXXX" );

    // Set parameters of elementary matrix
    elemMatr->setModel( _contact->getModel() );
    elemMatr->prepareCompute( option );

    // Prepare computing
    CalculPtr _calcul = std::make_unique< Calcul >( option );
    _calcul->setFiniteElementDescriptor( _fed );

    // Add input fields
    _calcul->addInputField( "PSTATCO", std::make_shared< FieldOnNodesReal >() );
    _calcul->addInputField( "PCOEFCO", std::make_shared< FieldOnNodesReal >() );
    _calcul->addInputField( "PDEPLPR", std::make_shared< FieldOnNodesReal >() );

    // Add output elementary
    _calcul->addOutputElementaryTerm( "PMATUUR", std::make_shared< ElementaryTermReal >() );

    // Compute
    _calcul->compute();

    elemMatr->addElementaryTerm( _calcul->getOutputElementaryTerm( "PMATUUR" ) );
    elemMatr->build();

    return elemMatr;
};