
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
#include "Utilities/Tools.h"

std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr >
ContactComputation::geometricGap( const ContactPairingPtr pairing ) const {

    // Prepare computing
    const std::string option( "GAP_GEOM" );
    CalculPtr calcul = std::make_unique< Calcul >( option );

    calcul->setFiniteElementDescriptor( pairing->getFiniteElementDescriptor() );

    // Add input fields
    calcul->addInputField( "PGEOMCR", pairing->getCoordinates() );

    // Add output elementary
    auto gap_elem = std::make_shared< ElementaryVectorReal >();
    gap_elem->prepareCompute( option );

    auto igap_elem = std::make_shared< ElementaryVectorReal >();
    igap_elem->prepareCompute( option );

    calcul->addOutputElementaryTerm( "PVECGAP", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVEIGAP", std::make_shared< ElementaryTermReal >() );

    // compute
    calcul->compute();

    // get elementary Term
    gap_elem->addElementaryTerm( calcul->getOutputElementaryTerm( "PVECGAP" ) );
    gap_elem->build();

    igap_elem->addElementaryTerm( calcul->getOutputElementaryTerm( "PVEIGAP" ) );
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
FieldOnCellsRealPtr ContactComputation::contactData( const ContactPairingPtr pairing,
                                                     const bool &initial_contact ) const {

    CALL_JEMARQ();

    auto fed = pairing->getFiniteElementDescriptor();

    // Field for intersection points and other thing ...
    auto data = std::make_shared< FieldOnCellsReal >( fed, "CHAR_MECA_CONT", "PCONFR" );

    ASTERINTEGER nbContPair = pairing->getNumberOfPairs();
    auto nbInter = concatenate( pairing->getNumberOfIntersectionPoints() );
    auto inter = concatenate( pairing->getSlaveIntersectionPoints() );
    AS_ASSERT( nbContPair == nbInter.size() );
    AS_ASSERT( 16 * nbContPair == inter.size() );

    auto pair2Zone = pairing->pairsToZones();

    auto grel = fed->getListOfGroupOfElements();
    grel->build();
    auto nbGrel = data->getNumberOfGroupOfElements();

    ASTERINTEGER nbPair = 0;

    // Loop on Grel
    for ( ASTERINTEGER iGrel = 0; iGrel < nbGrel; iGrel++ ) {
        auto nbElem = data->getNumberOfElements( iGrel );
        // size from catalogue
        AS_ASSERT( data->getSizeOfFieldOfElement( iGrel ) == 60 );
        auto liel = ( *grel )[iGrel + 1];
        for ( ASTERINTEGER iElem = 0; iElem < nbElem; iElem++ ) {
            auto iPair = -liel[iElem];

            if ( iPair <= nbContPair ) {
                auto iZone = pair2Zone[iPair - 1];
                auto zone = _contact->getContactZone( iZone );
                // Adress in field
                auto shift = data->getShifting( iGrel, iElem );

                // Number of intersection points
                ( *data )[shift + 0] = nbInter[iPair - 1];
                // Parametric coordinates
                for ( ASTERINTEGER i = 0; i < 16; i++ ) {
                    ( *data )[shift + 1 + i] = inter[16 * ( iPair - 1 ) + i];
                }

                /// Contact parameter
                auto cont = zone->getContactParameter();
                //  Value for ALGO_CONT
                ( *data )[shift + 23] = double( cont->getAlgorithm() );
                //  Value for TYPE_CONT
                ( *data )[shift + 24] = double( cont->getType() );
                //  Value for VARIANTE
                ( *data )[shift + 25] = double( cont->getVariant() );

                /// Friction parameter
                auto fric = zone->getFrictionParameter();
                //  Value for FROTTEMENT
                ( *data )[shift + 30] = fric->hasFriction();
                //  Value for ALGO_FROT
                ( *data )[shift + 31] = double( fric->getAlgorithm() );
                //  Value for TYPE_FROT
                ( *data )[shift + 32] = double( fric->getType() );
                // Value for threshold
                if ( fric->getType() == FrictionType::Tresca ) {
                    //  Value for TRESCA
                    ( *data )[shift + 34] = fric->getTresca();
                } else if ( fric->getType() == FrictionType::Coulomb ) {
                    //  Value for COULOMB
                    ( *data )[shift + 34] = fric->getCoulomb();
                }

                /// Other
                auto pair = zone->getPairingParameter();
                // Value for projection tolerancetolerance
                ( *data )[shift + 40] = 1.e-8;
                // Status to impose to contact
                if ( initial_contact ) {
                    ( *data )[shift + 41] = double( pair->getInitialState() );
                } else {
                    ( *data )[shift + 41] = double( InitialState::Interpenetrated );
                }

                nbPair++;
            }
        }
    }
    AS_ASSERT( nbPair == nbContPair );

    CALL_JEDEMA();

    return data;
};

/**
 * @brief Compute contact mortar matrix
 */
std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr >
ContactComputation::contactCoefficient() const {

    // Create FieldOnNodes
    auto ccont =
        std::make_shared< FieldOnNodesReal >( _contact->getFiniteElementDescriptor(), "ECCONT" );
    ccont->updateValuePointers();

    auto cfrot =
        std::make_shared< FieldOnNodesReal >( _contact->getFiniteElementDescriptor(), "ECFROT" );
    cfrot->updateValuePointers();

    auto dof2nodes = ccont->getNodesAndComponentsNumberFromDOF();
    std::map< ASTERINTEGER, ASTERINTEGER > nodes2dof;
    for ( ASTERINTEGER i_eq = 0; i_eq < dof2nodes.size(); i_eq++ ) {
        auto [node, cmp] = dof2nodes[i_eq];
        nodes2dof[node] = i_eq;
    }

    // Set values
    auto zones = _contact->getContactZones();

    for ( auto &zone : zones ) {
        auto coef_cont = zone->getContactParameter()->getCoefficient();

        auto coef_frot = zone->getFrictionParameter()->getCoefficient();
        auto nodes = zone->getSlaveNodes();
        for ( auto &node : nodes ) {
            ( *ccont )[nodes2dof[node]] = coef_cont;
            ( *cfrot )[nodes2dof[node]] = coef_frot;
        }
    }

    return std::make_pair( ccont, cfrot );
};
