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
      _nodeAndComponentsNumberFromDOF( getName() + ".DEEQ" ),
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

bool EquationNumbering::useLagrangeMultipliers() const {
    const std::string typeco( "NUME_EQUA" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXIS_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool EquationNumbering::useSingleLagrangeMultipliers() const {
    const std::string typeco( "NUME_EQUA" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "SIMP_LAGR" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

ASTERINTEGER
EquationNumbering::getNumberOfDofs( const bool local ) const {
    return _nodeAndComponentsNumberFromDOF->size() / 2;
};

void EquationNumbering::_buildAllComponentsNumber2Name() {
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
        auto name = trim( list_cmp[icmp - 1] );
        _componentsNumber2Name[icmp] = name;
        _componentsNumber2Name[-icmp] = "LAGR:" + name;
    }
};

VectorString EquationNumbering::getComponents() const {
    SetString ret;

    auto number2name = this->getComponentsNumber2Name();

    for ( auto &[num, name] : number2name ) {
        ret.insert( name );
    }

    return toVector( ret );
};

SetLong EquationNumbering::getComponentsNumber() const {
    auto ret = this->getNodesAndComponentsNumberFromDOF( true );

    SetLong cmpIds;
    for ( const auto &[nodeId, cmpId] : ret ) {
        cmpIds.insert( cmpId );
    }

    return cmpIds;
};

/**
 * @brief Maps between name of components and the nimber
 */
std::map< std::string, ASTERINTEGER > EquationNumbering::getComponentsName2Number() const {
    std::map< std::string, ASTERINTEGER > ret;

    auto number2name = this->getComponentsNumber2Name();
    for ( auto &[num, name] : number2name ) {
        ret[name] = num;
    }

    return ret;
};

std::map< ASTERINTEGER, std::string > EquationNumbering::getComponentsNumber2Name() const {
    std::map< ASTERINTEGER, std::string > ret;

    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsNumber2Name();

    auto cmpIds = this->getComponentsNumber();

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
    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDofs( true );

    VectorLong nodes( nb_eq );
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        // 0-based in c++
        nodes[i_eq] = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
    }

    return nodes;
}

VectorPairLong EquationNumbering::getNodesAndComponentsNumberFromDOF( const bool local ) const {
    const ASTERINTEGER nb_eq = this->getNumberOfDofs( true );

    VectorPairLong ret;
    ret.reserve( nb_eq );

    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        auto node_id = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
        auto cmp = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq + 1];
#ifdef ASTER_DEBUG_CXX
        AS_ASSERT( node_id >= 0 && node_id < getMesh()->getNumberOfNodes() );
#endif
        ret.push_back( std::make_pair( node_id, cmp ) );
    }

    return ret;
};

PairLong EquationNumbering::getNodeAndComponentNumberFromDOF( const ASTERINTEGER dof,
                                                              const bool local ) const {

    if ( dof < 0 or dof >= this->getNumberOfDofs( true ) ) {
        throw std::out_of_range( "Invalid node index: " + std::to_string( dof ) );
    }

    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    auto node_id = ( *_nodeAndComponentsNumberFromDOF )[2 * dof] - 1;
    auto cmp = ( *_nodeAndComponentsNumberFromDOF )[2 * dof + 1];

#ifdef ASTER_DEBUG_CXX
    AS_ASSERT( node_id >= 0 && node_id < getMesh()->getNumberOfNodes() );
#endif
    return std::make_pair( node_id, cmp );
};

std::vector< std::pair< ASTERINTEGER, std::string > >
EquationNumbering::getNodesAndComponentsFromDOF( const bool local ) const {
    auto nodesAndComponentsNumberFromDOF = this->getNodesAndComponentsNumberFromDOF( local );

    const ASTERINTEGER nb_eq = this->getNumberOfDofs( true );
    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsNumber2Name();

    std::vector< std::pair< ASTERINTEGER, std::string > > ret;
    ret.reserve( nodesAndComponentsNumberFromDOF.size() );

    for ( auto &[nodeId, cmpId] : nodesAndComponentsNumberFromDOF ) {
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
    auto [nodeId, cmpId] = this->getNodeAndComponentNumberFromDOF( dof, local );
    if ( _componentsNumber2Name.empty() )
        const_cast< EquationNumbering * >( this )->_buildAllComponentsNumber2Name();

#ifdef ASTER_DEBUG_CXX
    if ( _componentsNumber2Name.find( cmpId ) == _componentsNumber2Name.end() ) {
        std::cout << "Composante " << cmpId << " sans correspondance" << std::endl;
        raiseAsterError( "Erreur dans EquationNumbering" );
    }
#endif
    return std::make_pair( nodeId, _componentsNumber2Name.find( cmpId )->second );
};

std::map< PairLong, ASTERINTEGER >
EquationNumbering::getDOFsFromNodesAndComponentsNumber( const bool local ) const {
    auto descr = this->getNodesAndComponentsNumberFromDOF( local );

    std::map< PairLong, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
EquationNumbering::getDOFsFromNodesAndComponents( const bool local ) const {
    auto descr = this->getNodesAndComponentsFromDOF( local );

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

VectorLong EquationNumbering::getDOFs( const bool sameRank, const VectorString &list_cmp,
                                       const VectorLong &list_nodes ) const {
    const bool all_cmp = list_cmp.empty();
    const bool all_nodes = list_nodes.empty();
    const bool all_rank = !sameRank;
    const auto rank = getMPIRank();

    if ( !_mesh )
        raiseAsterError( "Mesh is empty" );

    const JeveuxVectorLong nodesRank = _mesh->getNodesRank();
    nodesRank->updateValuePointer();

    SetLong set_nodes = toSet( list_nodes );
    SetLong set_cmp;
    auto name2num = this->getComponentsName2Number();
    for ( auto &cmp : list_cmp ) {
        set_cmp.insert( name2num[trim( cmp )] );
    }

    auto descr = this->getNodesAndComponentsNumberFromDOF();
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

VectorLong EquationNumbering::getRowsAssociatedToPhysicalDofs( const bool local ) const {
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

VectorLong EquationNumbering::getRowsAssociatedToLagrangeMultipliers( const bool local ) const {
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

std::string EquationNumbering::getComponentAssociatedToRow( const ASTERINTEGER row,
                                                            const bool local ) const {
    auto [nodeId, cmpName] = this->getNodeAndComponentFromDOF( row );
    return cmpName;
};

ASTERINTEGER EquationNumbering::getNodeAssociatedToRow( const ASTERINTEGER row,
                                                        const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentNumberFromDOF( row );
    return nodeId;
};

bool EquationNumbering::isRowAssociatedToPhysical( const ASTERINTEGER row,
                                                   const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentNumberFromDOF( row );
    return cmpId > 0;
};

/**
 * @brief Mise a jour des pointeurs Jeveux
 * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
 */
void EquationNumbering::updateValuePointers() {
    _componentsOnNodes->build();
    _indexationVector->updateValuePointer();
    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    _lagrangianInformations->updateValuePointer();
};

bool EquationNumbering::operator==( EquationNumbering &toCompare ) {
    CALL_JEMARQ();
    bool ret = false;

    // TO FIX
    // if ( ( *_componentsOnNodes ) == ( *toCompare._componentsOnNodes ) ) {
    if ( ( *_indexationVector ) == ( *toCompare._indexationVector ) ) {
        if ( ( *_nodeAndComponentsNumberFromDOF ) ==
             ( *toCompare._nodeAndComponentsNumberFromDOF ) ) {
            ret = true;
        }
    }
    // }
    CALL_JEDEMA();

    return ret;
};
