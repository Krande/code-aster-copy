/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "Numbering/GlobalEquationNumbering.h"

#include "aster_fort_calcul.h"
#include "astercxx.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"
#include "Supervis/ResultNaming.h"

GlobalEquationNumbering::GlobalEquationNumbering( const std::string &baseName )
    : DataStructure( baseName, 19, "NUME_EQUA" ),
      _numberOfEquations( getName() + ".NEQU" ),
      _informations( getName() + ".REFN" ),
      _lagrangianInformations( getName() + ".DELG" ),
      _componentsOnNodes( getName() + ".PRNO" ),
      _namesOfGroupOfCells( getName() + ".LILI" ),
      _indexationVector( getName() + ".NUEQ" ),
      _nodeAndComponentsNumberFromDOF( getName() + ".DEEQ" ),
      _mesh( nullptr ),
      _model( nullptr ){};

void GlobalEquationNumbering::GlobalEquationNumbering::setModel( const ModelPtr &model ) {
    if ( _informations->exists() ) {
        _informations->updateValuePointer();
        const auto modelName = std::string( ( *_informations )[2].toString(), 0, 8 );
        if ( model && modelName != model->getName() )
            AS_ABORT( "Models are incompatible" );
    }
    _model = model;
    this->setMesh( _model->getMesh() );
};

void GlobalEquationNumbering::GlobalEquationNumbering::setMesh( const BaseMeshPtr &mesh ) {
    if ( _informations->exists() ) {
        _informations->updateValuePointer();
        const auto meshName = std::string( ( *_informations )[0].toString(), 0, 8 );
        if ( mesh && meshName != mesh->getName() )
            AS_ABORT( "Mesh are incompatible" );
    }
    _mesh = mesh;
};

std::string GlobalEquationNumbering::getPhysicalQuantity() const {
    _informations->updateValuePointer();
    JeveuxChar24 physicalQuantity = ( *_informations )[1];
    return physicalQuantity.rstrip();
};

ASTERINTEGER
GlobalEquationNumbering::getNumberOfDofs() const {
    return _nodeAndComponentsNumberFromDOF->size() / 2;
};

SetString GlobalEquationNumbering::getComponents() const {
    SetString ret;

    JeveuxVectorChar8 cmp_name( "&CMP_NAME" );
    JeveuxVectorLong cmp( "&CMP" );
    ASTERINTEGER ncmp;

    CALL_UTNCMP3( getName().c_str(), &ncmp, cmp->getName().c_str(), cmp_name->getName().c_str() );

    cmp_name->updateValuePointer();
    cmp->updateValuePointer();

    for ( int icmp = 0; icmp < ncmp; icmp++ ) {
        auto cmpId = ( *cmp )[icmp];
        std::string name;
        if ( cmpId > 0 ) {
            // Physical DoF
            name = trim( ( *cmp_name )[icmp] );
        } else if ( cmpId == 0 ) {
            // Lagrange multiplier associated to MPC
            name = "LAGR:MPC";
        } else {
            // Lagrange multiplier associated to Dirichlet BC
            name = "LAGR:" + trim( ( *cmp_name )[icmp] );
        }

        ret.insert( name );
    }

    return ret;
};

SetLong GlobalEquationNumbering::getComponentsNumber() const {

    JeveuxVectorChar8 cmp_name( "&CMP_NAME" );
    JeveuxVectorLong cmp( "&CMP" );
    ASTERINTEGER ncmp;

    CALL_UTNCMP3( getName().c_str(), &ncmp, cmp->getName().c_str(), cmp_name->getName().c_str() );

    return toSet( cmp->toVector() );
};

/**
 * @brief Maps between name of components and the nimber
 */
std::map< std::string, ASTERINTEGER > GlobalEquationNumbering::getComponentsName2Number() const {
    std::map< std::string, ASTERINTEGER > ret;

    JeveuxVectorChar8 cmp_name( "&CMP_NAME" );
    JeveuxVectorLong cmp( "&CMP" );
    ASTERINTEGER ncmp;

    CALL_UTNCMP3( getName().c_str(), &ncmp, cmp->getName().c_str(), cmp_name->getName().c_str() );

    cmp_name->updateValuePointer();
    cmp->updateValuePointer();

    for ( int icmp = 0; icmp < ncmp; icmp++ ) {
        auto cmpId = ( *cmp )[icmp];
        std::string name;
        if ( cmpId > 0 ) {
            // Physical DoF
            name = trim( ( *cmp_name )[icmp] );
        } else if ( cmpId == 0 ) {
            // Lagrange multiplier associated to MPC
            name = "LAGR:MPC";
        } else {
            // Lagrange multiplier associated to Dirichlet BC
            name = "LAGR:" + trim( ( *cmp_name )[icmp] );
        }

        ret[name] = cmpId;
    }

    return ret;
};

std::map< ASTERINTEGER, std::string > GlobalEquationNumbering::getComponentsNumber2Name() const {
    std::map< ASTERINTEGER, std::string > ret;

    auto name2number = this->getComponentsName2Number();
    for ( auto &[name, num] : name2number ) {
        ret[num] = name;
    }

    return ret;
};

VectorLong GlobalEquationNumbering::getNodesFromDOF() const {
    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDofs();

    VectorLong nodes( nb_eq );
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        // 0-based in c++
        nodes[i_eq] = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
    }

    return nodes;
}

VectorPairLong
GlobalEquationNumbering::getNodesAndComponentsNumberFromDOF( const bool local ) const {
    const ASTERINTEGER nb_eq = this->getNumberOfDofs();

    VectorPairLong ret;
    ret.reserve( nb_eq );

    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        auto node_id = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
        auto cmp = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq + 1];
        ret.push_back( std::make_pair( node_id, cmp ) );
    }

#ifdef ASTER_HAVE_MPI
    if ( !local ) {
        AS_ASSERT( _mesh );
        if ( _mesh->isParallel() ) {
            auto mapLG = _mesh->getLocalToGlobalMapping();
            mapLG->updateValuePointer();
            for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
                auto node_id = ret[i_eq].first;
                ret[i_eq].first = ( *mapLG )[node_id];
            }
        }
    }
#endif

    return ret;
};

PairLong GlobalEquationNumbering::getNodeAndComponentNumberFromDOF( const ASTERINTEGER dof,
                                                                    const bool local ) const {
    PairLong ret;
    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    const ASTERINTEGER nb_eq = this->getNumberOfDofs();
    for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
        auto node_id = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq] - 1;
        auto cmp = ( *_nodeAndComponentsNumberFromDOF )[2 * i_eq + 1];
        ret = std::make_pair( node_id, cmp );
    }

#ifdef ASTER_HAVE_MPI
    if ( !local ) {
        AS_ASSERT( _mesh );
        if ( _mesh->isParallel() ) {
            auto mapLG = _mesh->getLocalToGlobalMapping();
            mapLG->updateValuePointer();
            auto node_id = ret.first;
            ret.first = ( *mapLG )[node_id];
        }
    }
#endif

    return ret;
};

std::vector< std::pair< ASTERINTEGER, std::string > >
GlobalEquationNumbering::getNodesAndComponentsFromDOF( const bool local ) const {
    auto nodesAndComponentsNumberFromDOF = this->getNodesAndComponentsNumberFromDOF( local );

    const ASTERINTEGER nb_eq = this->getNumberOfDofs();
    auto num2name = this->getComponentsNumber2Name();

    std::vector< std::pair< ASTERINTEGER, std::string > > ret;
    ret.reserve( nb_eq );

    for ( auto &[nodeId, cmpId] : nodesAndComponentsNumberFromDOF ) {
        ret.push_back( std::make_pair( nodeId, num2name[cmpId] ) );
    }

    return ret;
};

std::pair< ASTERINTEGER, std::string >
GlobalEquationNumbering::getNodeAndComponentFromDOF( const ASTERINTEGER dof,
                                                     const bool local ) const {
    auto [nodeId, cmpId] = this->getNodeAndComponentNumberFromDOF( dof, local );

    auto num2name = this->getComponentsNumber2Name();

    return std::make_pair( nodeId, num2name[cmpId] );
};

std::map< PairLong, ASTERINTEGER >
GlobalEquationNumbering::getDOFsFromNodesAndComponentsNumber( const bool local ) const {
    auto descr = this->getNodesAndComponentsNumberFromDOF( local );

    std::map< PairLong, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER >
GlobalEquationNumbering::getDOFsFromNodesAndComponents( const bool local ) const {
    auto descr = this->getNodesAndComponentsFromDOF( local );

    std::map< std::pair< ASTERINTEGER, std::string >, ASTERINTEGER > ret;

    auto nbDof = descr.size();

    for ( ASTERINTEGER iDof = 0; iDof < nbDof; iDof++ ) {
        ret[descr[iDof]] = iDof;
    }

    return ret;
};

VectorLong GlobalEquationNumbering::getDOFs( const bool sameRank, const VectorString &list_cmp,
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

/**
 * @brief Mise a jour des pointeurs Jeveux
 * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
 */
void GlobalEquationNumbering::updateValuePointers() {
    _componentsOnNodes->build();
    _indexationVector->updateValuePointer();
    _nodeAndComponentsNumberFromDOF->updateValuePointer();
    _lagrangianInformations->updateValuePointer();
};

bool GlobalEquationNumbering::operator==( GlobalEquationNumbering &toCompare ) {
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
