
/**
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

#include "Contact/ContactComputation.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "DataFields/FieldOnNodes.h"
#include "Discretization/Calcul.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Messages/Messages.h"
#include "Numbering/DOFNumbering.h"
#include "Utilities/Tools.h"

FieldOnCellsRealPtr ContactComputation::contactData( const ContactPairingPtr contPairing,
                                                     const MaterialFieldPtr mater,
                                                     const bool &initial_contact ) const {
    CALL_JEMARQ();

    // Get finite element descriptor
    auto fed = contPairing->getFiniteElementDescriptor();

    // Field for intersection points and other thing ...
    auto data = FieldOnCellsPtrBuilder< ASTERDOUBLE >( fed, "CHAR_MECA_CONT", "PCONFR" );

    // Get pairing
    ASTERINTEGER nbContPair = contPairing->getNumberOfPairs();
    VectorLong nbInter = contPairing->getNumberOfIntersectionPoints();
    std::vector< VectorOfVectorsReal > inter =
        contPairing->getIntersectionPoints( CoordinatesSpace::Slave );
    VectorPairLong listPairs = contPairing->getListOfPairs();
    MapLong globPairToLocaPair = contPairing->globPairToLocaPair();

    AS_ASSERT( nbContPair == nbInter.size() );

    // Acces to list of cells
    const auto meshConnex = contPairing->getMesh()->getConnectivity();
    MapLong cellsToZones = contPairing->cellsToZones();
    auto grel = fed->getListOfGroupsOfElements();
    auto nbGrel = data->getNumberOfGroupOfElements();

    // Create local mapping (for Nitsche)
    auto mapping = []( const VectorLong &surf_nodes, const VectorLong &volu_nodes ) {
        VectorLong mapping;
        mapping.reserve( surf_nodes.size() );

        for ( auto &nodeId : surf_nodes ) {
            auto i = 1;
            for ( auto &nId : volu_nodes ) {
                if ( nId == nodeId ) {
                    break;
                }
                i++;
            }
            mapping.push_back( i );
        }

        return mapping;
    };

    // Get material (for Nitsche)
    auto listMaterial = mater->getVectorOfMaterial();

    // Loop on groups of elements
    ASTERINTEGER nbPair = 0;
    for ( ASTERINTEGER iGrel = 0; iGrel < nbGrel; iGrel++ ) {
        auto nbElem = data->getNumberOfElements( iGrel );
        AS_ASSERT( data->getSizeOfFieldOfElement( iGrel ) == 60 );
        auto liel = ( *grel )[iGrel + 1];
        liel->updateValuePointer();
        // Loop on elements
        for ( ASTERINTEGER iElem = 0; iElem < nbElem; iElem++ ) {
            // Get mesh cell index
            auto iPair = -( *liel )[iElem];

            // Current contact zone
            auto iZone = cellsToZones[iPair - 1];
            auto zone = _contact->getContactZone( iZone );

            // Adress in field
            auto shift = data->getShifting( iGrel, iElem );

            // Contact parameters
            auto cont = zone->getContactParameter();

            // Friction parameters
            auto fric = zone->getFrictionParameter();

            // Pairing parameters
            auto pair = zone->getPairingParameter();

            if ( iPair <= nbContPair ) {
                // Set number of intersection points
                ( *data )[shift + 0] = nbInter[iPair - 1];
                AS_ASSERT( nbInter[iPair - 1] <= 8 );

                // Set coordinates of slave intersection points
                for ( ASTERINTEGER iInter = 0; iInter < 8; iInter++ ) {
                    if ( iInter < nbInter[iPair - 1] ) {
                        auto iLocaPair = iPair - 1 - globPairToLocaPair[iZone];
                        ( *data )[shift + 1 + iInter] = inter[iZone][iLocaPair][2 * iInter];
                        ( *data )[shift + 9 + iInter] = inter[iZone][iLocaPair][2 * iInter + 1];
                    }
                }

                // For Nitsche
                if ( cont->getAlgorithm() == ContactAlgo::Nitsche ) {
                    auto [slavCellNume, mastCellNume] = listPairs[iPair - 1];
                    auto slav_surf_con = ( *meshConnex )[slavCellNume + 1]->toVector();
                    auto slavVoluNume = zone->getSlaveCellSurfToVolu( slavCellNume );
                    auto slav_volu_con = ( *meshConnex )[slavVoluNume + 1]->toVector();

                    auto mapLoc = mapping( slav_surf_con, slav_volu_con );
                    // Number of nodes
                    ( *data )[shift + 50] = double( mapLoc.size() );

                    // Mapping
                    auto i = 0;
                    for ( auto &nodeId : mapLoc ) {
                        ( *data )[shift + 51 + i++] = double( nodeId );
                    }

                    // Provisoire
                    AS_ASSERT( listMaterial.size() == 1 );
                    // Young modulus
                    ( *data )[shift + 45] = listMaterial[0]->getValueReal( "ELAS", "E" );
                    // Poisson ratio
                    ( *data )[shift + 46] = listMaterial[0]->getValueReal( "ELAS", "NU" );
                }

                // Value for projection tolerance
                ( *data )[shift + 40] = 1.e-8;

                // Status to impose to contact
                if ( initial_contact ) {
                    ( *data )[shift + 41] = double( pair->getInitialState() );
                } else {
                    ( *data )[shift + 41] = double( InitialState::Interpenetrated );
                }
            }

            //  Value for ALGO_CONT
            ( *data )[shift + 23] = double( cont->getAlgorithm() );
            //  Value for TYPE_CONT
            ( *data )[shift + 24] = double( cont->getType() );
            //  Value for VARIANTE
            ( *data )[shift + 25] = double( cont->getVariant() );
            //  Value for TYPE_MATR_TANG
            ( *data )[shift + 26] = double( cont->getJacobianType() );

            //  Value for FROTTEMENT
            ( *data )[shift + 30] = fric->hasFriction();
            //  Value for ALGO_FROT
            ( *data )[shift + 31] = double( fric->getAlgorithm() );
            //  Value for TYPE_FROT
            ( *data )[shift + 32] = double( fric->getType() );
            // Value for coefficient of friction
            if ( fric->getType() == FrictionType::Tresca ) {
                //  Value for TRESCA
                ( *data )[shift + 34] = fric->getTresca();
            } else if ( fric->getType() == FrictionType::Coulomb ) {
                //  Value for COULOMB
                ( *data )[shift + 34] = fric->getCoulomb();
            }

            nbPair++;
        }
    }

    CALL_JEDEMA();

    return data;
};

std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr >
ContactComputation::contactCoefficient() const {

    // Create FieldOnNodes
    auto fieldCont =
        std::make_shared< FieldOnNodesReal >( _contact->getFiniteElementDescriptor(), "ECCONT" );
    fieldCont->updateValuePointers();

    auto fieldFric =
        std::make_shared< FieldOnNodesReal >( _contact->getFiniteElementDescriptor(), "ECFROT" );
    fieldFric->updateValuePointers();

    auto dof2nodes = fieldCont->getDescription()->getNodeAndComponentIdFromDOF();
    MapLong nodes2dof;
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
            ( *fieldCont )[nodes2dof[node]] = coef_cont;
            ( *fieldFric )[nodes2dof[node]] = coef_frot;
        }
    }

    return std::make_pair( fieldCont, fieldFric );
};
