/**
 * @file CouplingPairing.cxx
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

#include "Coupling/CouplingPairing.h"

#include "aster_fort_ds.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Messages/Messages.h"
#include "Utilities/Tools.h"

namespace {
// Description of a virtual coupling cell

const std::map< std::pair< std::string, std::string >, std::string > cplCellNits {
    { { "MEDPTR3", "MECA_2D_HHO1_F" }, "CN_T3S3_HHO1" },
    { { "MEDPTR6", "MECA_2D_HHO1_F" }, "CP_T6S3_HHO1" },
    { { "MEDPQU4", "MECA_2D_HHO1_F" }, "CP_Q4S3_HHO2" },
    { { "MEDPQU8", "MECA_2D_HHO1_F" }, "CP_Q8S3_HHO2" },
    { { "MEDPQU9", "MECA_2D_HHO1_F" }, "CP_Q9S3_HHO3" },
};

const std::map< std::pair< std::string, std::string >, std::string > cplCellPena {
    { { "MEPLSE2", "MECA_2D_HHO1_F" }, "CP_S2S3_HHO1" },
    { { "MEPLSE3", "MECA_2D_HHO1_F" }, "CP_S3S3_HHO1" },
    { { "MEPLSE2", "MECA_2D_HHO2_F" }, "CP_S2S3_HHO2" },
    { { "MEPLSE3", "MECA_2D_HHO2_F" }, "CP_S3S3_HHO2" },
    { { "MEPLSE2", "MECA_2D_HHO3_F" }, "CP_S2S3_HHO3" },
    { { "MEPLSE3", "MECA_2D_HHO3_F" }, "CP_S3S3_HHO3" },
};
} // namespace

CouplingPairing::CouplingPairing( const std::string name, const ModelPtr model,
                                  const ASTERINTEGER verbosity )
    : DataStructure( name, 8, "PAIRING_SD" ),
      _model( model ),
      _meshPairing( std::make_shared< MeshPairing >( getName() + ".APMA" ) ),
      _verbosity( verbosity ),
      _algo( CouplingMethod::Undefined ),
      _coef_pena( 0. ) {
    if ( !_model )
        raiseAsterError( "Mesh is empty" );

    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );

    _meshPairing->setMesh( this->getMesh() );
    _meshPairing->setVerbosity( getVerbosity() );
    _meshPairing->setMethod( PairingMethod::Legacy );
};

void CouplingPairing::check() const {
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
    _meshPairing->checkNormals( _model );
}

void CouplingPairing::setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

ASTERBOOL CouplingPairing::compute() {

    _meshPairing->build();
    this->check();

    // Get distance ratio
    auto dist_pairing = -1.0;

    // Tolerance for pairing (zero value for geometrical operations)
    const ASTERDOUBLE pair_tole = 1e-8;

    // Tolerance for removing intersection cells
    const ASTERDOUBLE area_tole = 1e-2;

    // Pairing
    const auto returnValue = _meshPairing->compute( dist_pairing, pair_tole, area_tole );

    if ( !returnValue ) {
        return returnValue;
    }

    // Build FED
    this->buildFiniteElementDescriptor();

    return true;
}

ASTERINTEGER CouplingPairing::getCplCellType( const std::string &slavCellTypeName,
                                              const std::string &mastCellTypeName ) const {
    std::string cplTypeName;
    std::cout << slavCellTypeName << ", " << mastCellTypeName << std::endl;

    if ( slavCellTypeName.find( "HHO" ) != std::string::npos ) {
        UTMESS( "F", "COUPLAGE_2" );
    }

    const auto key = std::make_pair( slavCellTypeName, mastCellTypeName );

    if ( _algo == CouplingMethod::Nitsche ) {
        auto it = cplCellNits.find( key );
        if ( it != cplCellNits.end() ) {
            cplTypeName = it->second;
        }
    } else if ( _algo == CouplingMethod::Penalization ) {
        auto it = cplCellPena.find( key );
        if ( it != cplCellPena.end() ) {
            cplTypeName = it->second;
        }
    } else {
        AS_ABORT( "Not implemented" );
    };

    if ( cplTypeName.empty() ) {
        UTMESS( "F", "COUPLAGE_1" );
    }

    // Get index of type of coupling cell
    return _model->getFiniteElementDescriptor()->getElemTypeNume( cplTypeName );
}

void CouplingPairing::createVirtualElemForCoupling(
    MapLong &cplElemType, const JeveuxContiguousCollectionLong meshConnectivity,
    std::vector< VectorLong > &listCplElem, std::vector< VectorPairLong > &listCplType,
    SetLong &slaveNodePaired, SetLong &slaveCellPaired ) {

    // cplElemType: the number of elements for a given type
    // listCplType: list of coupling cells attached to pair (cellType, iPair)
    // listCplElem: list of coupling cells

    AS_ASSERT( _algo != CouplingMethod::Undefined );

    ASTERINTEGER iCplPair = 0;

    // Get pairing
    const auto listOfPairs = _meshPairing->getListOfPairs();
    const auto nbPairs = this->getNumberOfPairs();

    // Link between surface and volume
    MapLong surf2Volu;
    if ( _algo == CouplingMethod::Nitsche ) {
        if ( getMesh()->isParallel() ) {
            UTMESS( "F", "CONTACT1_5" );
        }
        surf2Volu = _meshPairing->getSlaveCellsSurfToVolu();
        if ( surf2Volu.size() == 0 ) {
            UTMESS( "F", "CONTACT1_3" );
        }
    }

    // Create vector of (virtual) coupling cells
    VectorPairLong listCplTypeZone;
    listCplTypeZone.reserve( nbPairs );

    const auto fed = _model->getFiniteElementDescriptor();

    const auto typeFE = fed->getFiniteElementType();
    typeFE->updateValuePointer();

    // Loop on pairs in zone
    for ( int iPair = 0; iPair < nbPairs; iPair++ ) {
        // Get pairing
        auto [slavCellNume, mastCellNume] = listOfPairs[iPair];

        // Get cell slave to construct (is volumic cell for Nitsche)
        auto slavCellUsedNume = slavCellNume;
        if ( _algo == CouplingMethod::Nitsche ) {
            slavCellUsedNume = surf2Volu[slavCellNume];
        }

        // Get slave and master cell type
        const auto slavCellTypeName = fed->getElemTypeName( ( *typeFE )[slavCellUsedNume] );
        const auto mastCellTypeName = fed->getElemTypeName( ( *typeFE )[mastCellNume] );

        // Get index of type of coupling cell
        const ASTERINTEGER typeElemNume =
            this->getCplCellType( slavCellTypeName, mastCellTypeName );
        AS_ASSERT( typeElemNume != -1 );

        // Add coupling element to list
        if ( cplElemType.count( typeElemNume ) == 0 ) {
            cplElemType[typeElemNume] = 0;
        }
        cplElemType[typeElemNume] += 1;

        // New virtual element
        iCplPair++;
        listCplTypeZone.push_back( std::make_pair( typeElemNume, iCplPair ) );

        // Get nodes
        auto slav_cell_con = ( *meshConnectivity )[slavCellUsedNume + 1];
        auto toAdd1 = slav_cell_con->toVector();
        slav_cell_con = JeveuxCollectionObject< ASTERINTEGER >();
        auto mast_cell_con = ( *meshConnectivity )[mastCellNume + 1];
        auto toAdd2 = mast_cell_con->toVector();

        // Coupling element on zone
        VectorLong cplElemZone;
        cplElemZone.reserve( toAdd1.size() + toAdd2.size() + 1 );

        // Copy slave nodes to coupling element
        cplElemZone.insert( cplElemZone.end(), toAdd1.begin(), toAdd1.end() );

        // Add slave nodes to list of paired nodes
        slaveNodePaired.insert( toAdd1.begin(), toAdd1.end() );

        // Add slave cell to list of paired cells
        slaveCellPaired.insert( slavCellNume );

        // Copy master nodes to coupling element
        cplElemZone.insert( cplElemZone.end(), toAdd2.begin(), toAdd2.end() );

        // Add type of coupling element
        cplElemZone.push_back( typeElemNume );

        // Add coupling element to all coupling elements
        listCplElem.push_back( cplElemZone );
    }
    // Add coupling elements of zone
    listCplType.push_back( listCplTypeZone );

    AS_ASSERT( iCplPair == nbPairs )
}

void CouplingPairing::buildFiniteElementDescriptor() {

    CALL_JEMARQ();

    // Mesh
    auto mesh = getMesh();
    const auto meshConnectivity = mesh->getConnectivity();

    // Get pairing parameters
    const ASTERINTEGER nbCplPairTot = this->getNumberOfPairs();

    // Create objets for nodes and cells
    VectorOfVectorsLong listCplElem;
    listCplElem.reserve( nbCplPairTot );
    std::vector< VectorPairLong > listCplType;
    listCplType.reserve( 2 );

    // Objects
    SetLong slaveNodePaired, slaveCellPaired;

    // Object for number of cells for each type of coupling cell
    MapLong cplElemType;

    // Create virtual elements for coupling
    createVirtualElemForCoupling( cplElemType, meshConnectivity, listCplElem, listCplType,
                                  slaveNodePaired, slaveCellPaired );

    // Create finite element descriptor for virtual coupling elements
    _fed = std::make_shared< FiniteElementDescriptor >( mesh );

    // Create list of virtual nodes (none !)
    _fed->setNumberOfVirtualNodes( 0 );

    // Create list of virtual elements (NEMA object)
    auto cplResFEDNema = _fed->getVirtualCellsDescriptor();
    cplResFEDNema->allocate( listCplElem );

    // Number of groups of elements and length of FED for coupling element
    ASTERINTEGER nbGrel = 0, lielLont = 0;
    for ( auto &[type, size] : cplElemType ) {
        lielLont += size;
        nbGrel += 1;
    }
    lielLont += nbGrel;

    // Create list of elements (LIEL object)
    auto cplResFEDLiel = _fed->getListOfGroupsOfElements();
    cplResFEDLiel->allocate( nbGrel, lielLont, Variable );

    // Add virtual elements for each GREL
    for ( auto &[type, size] : cplElemType ) {
        VectorLong virtualCells;
        virtualCells.reserve( size );

        for ( auto &listCplTypeZone : listCplType ) {
            for ( auto &[typeElemNume, iCplPair] : listCplTypeZone ) {
                if ( typeElemNume == type ) {
                    // Virtual cells = index of cells is negative
                    virtualCells.push_back( -iCplPair );
                }
            }
        }
        virtualCells.push_back( type );
        cplResFEDLiel->push_back( virtualCells );
    }

    // Get parameters from standard model
    auto paramToCopy = _model->getFiniteElementDescriptor()->getParameters();
    paramToCopy->updateValuePointer();
    auto docu = paramToCopy->getInformationParameter();

    // Create LGRF object
    auto parameters = _fed->getParameters();
    parameters->allocate( 4 );
    ( *parameters )[0] = mesh->getName();
    ( *parameters )[1] = _model->getName();
    ( *parameters )[2] = ( *paramToCopy )[2];
    parameters->setInformationParameter( docu );

    // Adapt FED
    CALLO_ADALIG_WRAP( _fed->getName() );
    bool l_calc_rigi = false;
    CALLO_INITEL( _fed->getName(), (ASTERLOGICAL *)&l_calc_rigi );

    // Final building
    _fed->build();

    CALL_JEDEMA();
};

FieldOnCellsRealPtr CouplingPairing::getPairingField() const {
    CALL_JEMARQ();

    // Field for intersection points and other thing ...
    auto data = FieldOnCellsPtrBuilder< ASTERDOUBLE >( _fed, "CHAR_MECA_CPL", "PPAIRR" );

    // Get pairing
    const ASTERINTEGER nbPairs = this->getNumberOfPairs();
    const VectorLong nbInter = _meshPairing->getNumberOfIntersectionPoints();
    AS_ASSERT( nbPairs == nbInter.size() );

    VectorPairLong listPairs = _meshPairing->getListOfPairs();

    // Acces to list of cells
    const auto meshConnex = this->getMesh()->getConnectivity();
    auto grel = _fed->getListOfGroupsOfElements();
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

    // Loop on groups of elements
    ASTERINTEGER nbPair = 0;
    for ( ASTERINTEGER iGrel = 0; iGrel < nbGrel; iGrel++ ) {
        auto nbElem = data->getNumberOfElements( iGrel );
        AS_ASSERT( data->getSizeOfFieldOfElement( iGrel ) == 41 );
        auto liel = ( *grel )[iGrel + 1];
        liel->updateValuePointer();
        // Loop on elements
        for ( ASTERINTEGER iElem = 0; iElem < nbElem; iElem++ ) {
            // Get mesh cell index
            auto iPair = -( *liel )[iElem];

            // Adress in field
            auto shift = data->getShifting( iGrel, iElem );

            if ( iPair <= nbPairs ) {
                const auto iLocaPair = iPair - 1;

                // Value for projection tolerance
                ( *data )[shift + 0] = 1.e-8;
                // Set number of intersection points
                ( *data )[shift + 1] = nbInter[iLocaPair];
                AS_ASSERT( nbInter[iLocaPair] <= 8 );

                // Set coordinates of slave intersection points
                const auto inter =
                    _meshPairing->getIntersectionPoints( iLocaPair, CoordinatesSpace::Slave );

                for ( ASTERINTEGER iInter = 0; iInter < nbInter[iLocaPair]; iInter++ ) {
                    ( *data )[shift + 2 + iInter] = inter[iInter][0];
                    ( *data )[shift + 10 + iInter] = inter[iInter][1];
                }

                // Coefficient of penalization
                ( *data )[shift + 29] = _coef_pena;

                // For Nitsche
                if ( _algo == CouplingMethod::Nitsche ) {
                    auto [slavCellNume, mastCellNume] = listPairs[iLocaPair];
                    auto slav_surf_con = ( *meshConnex )[slavCellNume + 1]->toVector();
                    auto slavVoluNume = _meshPairing->getSlaveCellSurfToVolu( slavCellNume );
                    auto slav_volu_con = ( *meshConnex )[slavVoluNume + 1]->toVector();

                    auto mapLoc = mapping( slav_surf_con, slav_volu_con );
                    // Number of nodes
                    ( *data )[shift + 30] = double( mapLoc.size() );

                    // Mapping
                    auto i = 0;
                    for ( auto &nodeId : mapLoc ) {
                        ( *data )[shift + 31 + i++] = double( nodeId );
                    }
                }
            }

            nbPair++;
        }
    }

    CALL_JEDEMA();

    return data;
};
