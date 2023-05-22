/**
 * @file EquationNumbering.cxx
 * @brief Implementation de EquationNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "Numbering/EquationNumbering.h"

#include "aster_fort_calcul.h"

#include "Modeling/PhysicalQuantityManager.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/Tools.h"

EquationNumbering::EquationNumbering( const std::string &baseName )
    : BaseEquationNumbering( baseName ),
      _numberOfEquations( getName() + ".NEQU" ),
      _informations( getName() + ".REFN" ),
      _lagrangianInformations( getName() + ".DELG" ),
      _componentsOnNodes( getName() + ".PRNO" ),
      _namesOfGroupOfCells( getName() + ".LILI" ),
      _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsIdFromDOF( getName() + ".DEEQ" ),
      _mesh( nullptr ),
      _model( nullptr ) {};

EquationNumbering::EquationNumbering() : EquationNumbering( DataStructureNaming::getNewName() ) {};

bool EquationNumbering::exists() const {
    return _informations.exists() && _componentsOnNodes.exists() && _indexationVector.exists();
};

void EquationNumbering::EquationNumbering::setModel( const ModelPtr &model ) {
    if ( model && exists() ) {
        _informations->updateValuePointer();
        const auto modelName = std::string( ( *_informations )[2].toString(), 0, 8 );
        if ( model && modelName != model->getName() ) {
            AS_ABORT( "Models are incompatible" );
        }
        _model = model;
        this->setMesh( _model->getMesh() );
    }
};

void EquationNumbering::EquationNumbering::setMesh( const BaseMeshPtr &mesh ) {
    if ( mesh && exists() ) {
        _informations->updateValuePointer();
        const auto meshName = std::string( ( *_informations )[0].toString(), 0, 8 );
        if ( mesh && meshName != mesh->getName() ) {
            raiseAsterError( "Mesh are incompatible: " + mesh->getName() + " vs " + meshName );
        }
        _mesh = mesh;
    }
};

std::string EquationNumbering::getPhysicalQuantity() const {
    _informations->updateValuePointer();
    JeveuxChar24 physicalQuantity = ( *_informations )[1];
    return physicalQuantity.rstrip();
};

bool EquationNumbering::useLagrangeDOF() const {
    const std::string typeco( "NUME_EQUA" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXIS_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = strip( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool EquationNumbering::useSingleLagrangeDOF() const {
    const std::string typeco( "NUME_EQUA" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "SIMP_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = strip( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

ASTERINTEGER
EquationNumbering::getNumberOfDOFs( const bool local ) const {
    return _nodeAndComponentsIdFromDOF->size() / 2;
};

void EquationNumbering::_buildAllComponentsId2Name() {
    const std::string typeco( "NUME_EQUA" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "F" );
    const std::string questi( "NUM_GD" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

    auto list_cmp = PhysicalQuantityManager::getComponentNames( repi );
    int nb_cmp = list_cmp.size();

    _componentsNumber2Name[0] = "LAGR:MPC";
    for ( int icmp = 1; icmp <= nb_cmp; icmp++ ) {
        auto name = strip( list_cmp[icmp - 1] );
        _componentsNumber2Name[icmp] = name;
        _componentsNumber2Name[-icmp] = "LAGR:" + name;
    }
};

VectorString EquationNumbering::getComponents() const {
    SetString ret;

    auto number2name = this->getComponentsIdToName();

    for ( auto &[num, name] : number2name ) {
        ret.insert( name );
    }

    return toVector( ret );
};

SetLong EquationNumbering::getComponentsId() const {
    auto ret = this->getNodeAndComponentIdFromDOF( true );

    SetLong cmpIds;
    for ( const auto &[nodeId, cmpId] : ret ) {
        cmpIds.insert( cmpId );
    }

    return cmpIds;
};

/**
 * @brief Maps between name of components and the nimber
 */
std::map< std::string, ASTERINTEGER > EquationNumbering::getComponentsNameToId() const {
    std::map< std::string, ASTERINTEGER > ret;

    auto number2name = this->getComponentsIdToName();
    for ( auto &[num, name] : number2name ) {
        ret[name] = num;
    }

    return ret;
};

std::map< ASTERINTEGER, std::string > EquationNumbering::getComponentsIdToName() const {
    std::map< ASTERINTEGER, std::string > ret;

    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsId2Name();

    auto cmpIds = this->getComponentsId();

    for ( auto &cmpId : cmpIds ) {
        ret[cmpId] = _componentsNumber2Name.find( cmpId )->second;

#ifdef ASTER_DEBUG_CXX
        if ( _componentsNumber2Name.find( cmpId ) == _componentsNumber2Name.end() ) {
            std::cout << "Composante " << cmpId << " sans correspondance" << std::endl;
            raiseAsterError( "Erreur dans EquationNumbering" );
        }
#endif
    }

    return ret;
};

VectorLong EquationNumbering::getNodesFromDOF() const {
    _nodeAndComponentsIdFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDOFs( true );

    VectorLong nodes( nb_eq );
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        // 0-based in c++
        nodes[i_eq] = ( *_nodeAndComponentsIdFromDOF )[2 * i_eq] - 1;
    }

    return nodes;
}

VectorPairLong EquationNumbering::getNodeAndComponentIdFromDOF( const bool local ) const {
    const ASTERINTEGER nb_eq = this->getNumberOfDOFs( true );

    VectorPairLong ret;
    ret.reserve( nb_eq );

    _nodeAndComponentsIdFromDOF->updateValuePointer();
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        auto node_id = ( *_nodeAndComponentsIdFromDOF )[2 * i_eq] - 1;
        auto cmp = ( *_nodeAndComponentsIdFromDOF )[2 * i_eq + 1];
#ifdef ASTER_DEBUG_CXX
        AS_ASSERT( node_id < 0 || node_id < getMesh()->getNumberOfNodes() );
#endif
        ret.push_back( std::make_pair( node_id, cmp ) );
    }

    return ret;
};

PairLong EquationNumbering::getNodeAndComponentIdFromDOF( const ASTERINTEGER dof,
                                                          const bool local ) const {

    if ( dof < 0 or dof >= this->getNumberOfDOFs( true ) ) {
        throw std::out_of_range( "Invalid node index: " + std::to_string( dof ) );
    }

    _nodeAndComponentsIdFromDOF->updateValuePointer();
    auto node_id = ( *_nodeAndComponentsIdFromDOF )[2 * dof] - 1;
    auto cmp = ( *_nodeAndComponentsIdFromDOF )[2 * dof + 1];

#ifdef ASTER_DEBUG_CXX
    AS_ASSERT( node_id < 0 || node_id < getMesh()->getNumberOfNodes() );
#endif
    return std::make_pair( node_id, cmp );
};

std::vector< std::pair< ASTERINTEGER, std::string > >
EquationNumbering::getNodeAndComponentFromDOF( const bool local ) const {
    auto nodesAndComponentsIdFromDOF = this->getNodeAndComponentIdFromDOF( local );

    const ASTERINTEGER nb_eq = this->getNumberOfDOFs( true );
    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsId2Name();

    std::vector< std::pair< ASTERINTEGER, std::string > > ret;
    ret.reserve( nodesAndComponentsIdFromDOF.size() );

    for ( auto &[nodeId, cmpId] : nodesAndComponentsIdFromDOF ) {
        ret.push_back( std::make_pair( nodeId, _componentsNumber2Name.find( cmpId )->second ) );

#ifdef ASTER_DEBUG_CXX
        if ( _componentsNumber2Name.find( cmpId ) == _componentsNumber2Name.end() ) {
            std::cout << "Composante " << cmpId << " sans correspondance" << std::endl;
            raiseAsterError( "Erreur dans EquationNumbering" );
        }
#endif
    }

    return ret;
};

std::pair< ASTERINTEGER, std::string >
EquationNumbering::getNodeAndComponentFromDOF( const ASTERINTEGER dof, const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentIdFromDOF( dof, local );
    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsId2Name();

#ifdef ASTER_DEBUG_CXX
    if ( _componentsNumber2Name.find( cmpId ) == _componentsNumber2Name.end() ) {
        std::cout << "Composante " << cmpId << " sans correspondance" << std::endl;
        raiseAsterError( "Erreur dans EquationNumbering" );
    }
#endif
    return std::make_pair( nodeId, _componentsNumber2Name.find( cmpId )->second );
};

std::map< PairLong, ASTERINTEGER >
EquationNumbering::getDOFFromNodeAndComponentId( const bool local ) const {
    auto descr = this->getNodeAndComponentIdFromDOF( local );

    std::map< PairLong, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
EquationNumbering::getDOFFromNodeAndComponent( const bool local ) const {
    auto descr = this->getNodeAndComponentFromDOF( local );

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

ASTERINTEGER EquationNumbering::getDOFFromNodeAndComponent( const ASTERINTEGER &node,
                                                            const std::string &comp,
                                                            const bool local ) const {
    auto dofs = this->getDOFFromNodeAndComponent( local );

    for ( auto it = dofs.cbegin(); it != dofs.cend(); ++it ) {
        auto [nodeId, cmp] = it->first;
        if ( node == nodeId && comp == cmp ) {
            return it->second;
        }
    }

    throw std::runtime_error( "DOF nout found" );

    return -1;
};

VectorLong EquationNumbering::getDOFs( const bool sameRank, const VectorString &list_cmp,
                                       const VectorString &list_grpno ) const {
    const bool all_cmp = list_cmp.empty();
    const bool all_nodes = list_grpno.empty();
    const bool all_rank = !sameRank;
    const auto rank = getMPIRank();

    if ( !_mesh )
        raiseAsterError( "Mesh is empty" );

    const JeveuxVectorLong nodesRank = _mesh->getNodesRank();
    nodesRank->updateValuePointer();

    SetLong set_nodes = toSet( getMesh()->getNodes( list_grpno ) );
    SetLong set_cmp;
    auto name2num = this->getComponentsNameToId();
    for ( auto &cmp : list_cmp ) {
        set_cmp.insert( name2num[strip( cmp )] );
    }

    auto descr = this->getNodeAndComponentIdFromDOF();
    auto nbDof = descr.size();

    VectorLong dofUsed;
    dofUsed.reserve( nbDof );

    for ( auto dof = 0; dof < nbDof; ++dof ) {
        const auto node_id = std::abs( descr[dof].first );
        const auto cmp_id = descr[dof].second;
        const bool l_keep_cmp = all_cmp || ( set_cmp.count( cmp_id ) > 0 );
        const bool l_keep_rank = all_rank || ( ( *nodesRank )[node_id] == rank );
        const bool l_keep_node = all_nodes || ( set_nodes.count( node_id ) > 0 );
        if ( l_keep_node && l_keep_cmp && l_keep_rank ) {
            dofUsed.push_back( dof );
        }
    }

    return dofUsed;
};

VectorLong EquationNumbering::getPhysicalDOFs( const bool local ) const {
    auto lagrInfo = this->getLagrangianInformations();
    lagrInfo->updateValuePointer();
    ASTERINTEGER size = lagrInfo->size();
    VectorLong physicalRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *lagrInfo )[i];
        if ( physicalIndicator == 0 )
            physicalRows.push_back( i );
    }
    return physicalRows;
};

VectorLong EquationNumbering::getLagrangeDOFs( const bool local ) const {
    auto lagrInfo = this->getLagrangianInformations();
    lagrInfo->updateValuePointer();
    ASTERINTEGER size = lagrInfo->size();
    VectorLong lagrangeRows;
    ASTERINTEGER physicalIndicator;
    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *lagrInfo )[i];
        if ( physicalIndicator != 0 )
            lagrangeRows.push_back( i );
    }
    return lagrangeRows;
};

std::string EquationNumbering::getComponentFromDOF( const ASTERINTEGER dof,
                                                    const bool local ) const {
    auto [nodeId, cmpName] = this->getNodeAndComponentFromDOF( dof );
    return cmpName;
};

ASTERINTEGER EquationNumbering::getNodeFromDOF( const ASTERINTEGER dof, const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentIdFromDOF( dof );
    return nodeId;
};

bool EquationNumbering::isPhysicalDOF( const ASTERINTEGER dof, const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentIdFromDOF( dof );
    return cmpId > 0;
};

std::pair< std::pair< VectorLong, VectorString >, VectorLong >
EquationNumbering::getDOFsWithDescription( const std::string cmp,
                                           const VectorString groupNames ) const {

    VectorLong v_nodes;
    VectorString cmps;
    VectorLong dofs;

    std::set< ASTERINTEGER > nodes;
    if ( groupNames.size() == 0 ) {
        auto group = _mesh->getNodes();
        std::copy( group.begin(), group.end(), std::inserter( nodes, nodes.end() ) );
    } else {
        for ( auto groupName : groupNames ) {
            auto group = _mesh->getNodes( groupName );
            std::copy( group.begin(), group.end(), std::inserter( nodes, nodes.end() ) );
        }
    }

    auto idToName = getComponentsIdToName();

    ASTERINTEGER icmp, ncmp;
    bool all_cmp = strip( cmp ) == "";
    if ( all_cmp ) {
        ncmp = getComponents().size();
        cmps.reserve( ncmp * nodes.size() );
    } else {
        ncmp = 1;
        icmp = getComponentsNameToId()[cmp];
    }
    v_nodes.reserve( ncmp * nodes.size() );

    const auto descr = getNodeAndComponentIdFromDOF();

    for ( auto dof = 0; dof < descr.size(); ++dof ) {
        if ( all_cmp and descr[dof].second > 0 ) {
            if ( nodes.find( descr[dof].first ) != nodes.end() ) {
                v_nodes.push_back( descr[dof].first );
                cmps.push_back( idToName[descr[dof].second] );
                dofs.push_back( dof );
            }
        } else if ( icmp == descr[dof].second ) {
            if ( nodes.find( descr[dof].first ) != nodes.end() ) {
                v_nodes.push_back( descr[dof].first );
                dofs.push_back( dof );
            }
        }
    }

    return std::make_pair( std::make_pair( v_nodes, cmps ), dofs );
};

/**
 * @brief Mise a jour des pointeurs Jeveux
 * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
 */
void EquationNumbering::updateValuePointers() {
    _componentsOnNodes->build();
    _indexationVector->updateValuePointer();
    _nodeAndComponentsIdFromDOF->updateValuePointer();
    _lagrangianInformations->updateValuePointer();
};

bool EquationNumbering::operator==( EquationNumbering &toCompare ) {
    CALL_JEMARQ();
    bool ret = false;

    // TO FIX
    // if ( ( *_componentsOnNodes ) == ( *toCompare._componentsOnNodes ) ) {
    if ( ( *_indexationVector ) == ( *toCompare._indexationVector ) ) {
        if ( ( *_nodeAndComponentsIdFromDOF ) == ( *toCompare._nodeAndComponentsIdFromDOF ) ) {
            ret = true;
        }
    }
    // }
    CALL_JEDEMA();

    return ret;
};
