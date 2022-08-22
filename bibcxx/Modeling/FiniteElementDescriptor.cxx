/**
 * @file FiniteElementDescriptor.cxx
 * @brief Implementation de FiniteElementDescriptor
 * @author Nicolas Sellenet
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

#include "Modeling/FiniteElementDescriptor.h"

#include "aster_fort_ds.h"

#include "Meshes/BaseMesh.h"
#include "Meshes/ConnectionMesh.h"
#include "Modeling/PhysicalQuantityManager.h"
#include "Modeling/PhysicsAndModelings.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

FiniteElementDescriptor::FiniteElementDescriptor( const std::string &name, const BaseMeshPtr mesh )
    : DataStructure( name, 19, "LIGREL" ),
      _numberOfDelayedNumberedConstraintNodes( getName() + ".NBNO" ),
      _parameters( getName() + ".LGRF" ),
      _dofDescriptor( getName() + ".PRNM" ),
      _listOfGroupOfCells( getName() + ".LIEL" ),
      _groupsOfCellsNumberByElement( getName() + ".REPE" ),
      _delayedNumberedConstraintElementsDescriptor( getName() + ".NEMA" ),
      _dofOfDelayedNumberedConstraintNodes( getName() + ".PRNS" ),
      _virtualNodesNumbering( getName() + ".LGNS" ),
      _superElementsDescriptor( getName() + ".SSSA" ),
      _nameOfNeighborhoodStructure( getName() + ".NVGE" ),
      _mesh( mesh ),
      _explorer( FiniteElementDescriptor::ConnectivityVirtualCellsExplorer(
          _delayedNumberedConstraintElementsDescriptor ) ),
      _explorer2(
          FiniteElementDescriptor::ConnectivityVirtualCellsExplorer( _listOfGroupOfCells ) ){};

FiniteElementDescriptor::FiniteElementDescriptor( const BaseMeshPtr mesh )
    : FiniteElementDescriptor( DataStructureNaming::getNewName(), mesh ){};

FiniteElementDescriptor::FiniteElementDescriptor( const ModelPtr model )
    : FiniteElementDescriptor( model->getMesh() ){};

FiniteElementDescriptor::FiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc,
                                                  const VectorString &groupOfCells )
    : FiniteElementDescriptor( FEDesc->getMesh() ) {

    VectorLong commonCells;

    for ( auto &group : groupOfCells ) {
        auto cells = _mesh->getCells( group );
        auto it = commonCells.end();
        commonCells.insert( it, cells.begin(), cells.end() );
    }

    VectorLong listOfCells = unique( commonCells );

    std::string base( "G" );
    ASTERINTEGER nbCells = listOfCells.size();
    for ( auto &cell : listOfCells )
        cell += 1;
    CALL_EXLIM2( listOfCells.data(), &nbCells, FEDesc->getName(), base, getName() );
};

FiniteElementDescriptor::FiniteElementDescriptor( const ModelPtr model,
                                                  const VectorString &groupOfCells )
    : FiniteElementDescriptor( model->getFiniteElementDescriptor(), groupOfCells ) {
    setModel( model );
}

const FiniteElementDescriptor::ConnectivityVirtualCellsExplorer &
FiniteElementDescriptor::getVirtualCellsExplorer() const {
    _delayedNumberedConstraintElementsDescriptor->build();
    return _explorer;
};

const JeveuxVectorLong &FiniteElementDescriptor::getVirtualNodesComponentDescriptor() const {
    return _dofOfDelayedNumberedConstraintNodes;
};

const JeveuxVectorLong &FiniteElementDescriptor::getVirtualNodesNumbering() const {
    _virtualNodesNumbering->updateValuePointer();
    return _virtualNodesNumbering;
};

const FiniteElementDescriptor::ConnectivityVirtualCellsExplorer &
FiniteElementDescriptor::getListOfGroupOfElementsExplorer() const {
    _listOfGroupOfCells->build();
    return _explorer2;
};

const JeveuxCollectionLong &FiniteElementDescriptor::getListOfGroupOfElements() const {
    _listOfGroupOfCells->build();
    return _listOfGroupOfCells;
};

const JeveuxCollectionLong &FiniteElementDescriptor::getVirtualCellsDescriptor() const {
    _delayedNumberedConstraintElementsDescriptor->build();
    return _delayedNumberedConstraintElementsDescriptor;
};

ASTERINTEGER FiniteElementDescriptor::getNumberOfVirtualNodes() const {
    _numberOfDelayedNumberedConstraintNodes->updateValuePointer();
    return ( *_numberOfDelayedNumberedConstraintNodes )[0];
};

void FiniteElementDescriptor::setNumberOfVirtualNodes( const ASTERINTEGER nbNodes ) {

    if ( !_numberOfDelayedNumberedConstraintNodes->exists() ) {
        _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    } else {
        _numberOfDelayedNumberedConstraintNodes->updateValuePointer();
    }
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = nbNodes;
};

JeveuxVectorLong FiniteElementDescriptor::getNumberOfVirtualNodesDescriptor() {
    return _numberOfDelayedNumberedConstraintNodes;
};

JeveuxVectorChar8 FiniteElementDescriptor::getParameters() const { return _parameters; };

const JeveuxVectorLong &FiniteElementDescriptor::getPhysicalNodesComponentDescriptor() const {
    _dofDescriptor->updateValuePointer();
    return _dofDescriptor;
};

const JeveuxVectorLong &FiniteElementDescriptor::getListOfGroupOfElementsbyElement() const {
    _groupsOfCellsNumberByElement->updateValuePointer();
    return _groupsOfCellsNumberByElement;
};

const BaseMeshPtr FiniteElementDescriptor::getMesh() const { return _mesh; };

void FiniteElementDescriptor::setMesh( const BaseMeshPtr &currentMesh ) { _mesh = currentMesh; };

int FiniteElementDescriptor::getPhysics( void ) const {
    const std::string docu = trim( _parameters->getInformationParameter() );

    if ( docu == "MECA" )
        return Physics::Mechanics;
    else if ( docu == "THER" )
        return Physics::Thermal;
    else if ( docu == "ACOU" )
        return Physics::Acoustic;
    else
        throw std::runtime_error( "Unknown physics" );

    return -1;
};

void FiniteElementDescriptor::setModel( const ModelPtr model ) { _model = model; };

ModelPtr FiniteElementDescriptor::getModel() { return _model.lock(); }

ASTERINTEGER FiniteElementDescriptor::numberOfSuperElement() {
    const std::string typeco( "LIGREL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "NB_SS_ACTI" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

    return repi;
};

bool FiniteElementDescriptor::existsFiniteElement() {
    const std::string typeco( "LIGREL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_ELEM" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};


bool FiniteElementDescriptor::existsSuperElement() { return ( numberOfSuperElement() > 0 ); }

bool FiniteElementDescriptor::exists() const {
    if ( _parameters->exists() && _numberOfDelayedNumberedConstraintNodes->exists() )
        return true;

    return false;
};

#ifdef ASTER_HAVE_MPI
void FiniteElementDescriptor::transferDofDescriptorFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    const JeveuxVectorLong &otherDofDescriptor = other->getPhysicalNodesComponentDescriptor();

    const int rank = getMPIRank();
    const int size = getMPISize();
    int nbNodes = connectionMesh->getNumberOfNodes();
    int nec = otherDofDescriptor->size() / other->getMesh()->getNumberOfNodes();

    const JeveuxVectorLong &localNumbering = connectionMesh->getNodesLocalNumbering();
    const JeveuxVectorLong &owner = connectionMesh->getNodesOwner();

    int nbNodesLoc = 0;
    for ( int i = 0; i < nbNodes; ++i ) {
        if ( ( *owner )[i] == rank )
            nbNodesLoc++;
    }

    VectorLong buffer( nec * nbNodesLoc, 0 );
    nbNodesLoc = 0;
    for ( int i = 0; i < nbNodes; ++i ) {
        int nodeNum = ( *localNumbering )[i] - 1;
        if ( ( *owner )[i] == rank ) {
            for ( int j = 0; j < nec; ++j )
                buffer[nbNodesLoc * nec + j] = ( *otherDofDescriptor )[nodeNum * nec + j];

            nbNodesLoc++;
        }
    }

    std::vector< VectorLong > gathered;
    AsterMPI::all_gather( buffer, gathered );
    buffer.clear();

    _dofDescriptor->allocate( nbNodes * nec );
    VectorLong nbNodesProc( size, 0 );
    for ( int i = 0; i < nbNodes; ++i ) {
        auto rowner = ( *owner )[i];
        for ( int j = 0; j < nec; ++j )
            ( *_dofDescriptor )[i * nec + j] = gathered[rowner][nbNodesProc[rowner] * nec + j];
        nbNodesProc[rowner]++;
    }
};

void FiniteElementDescriptor::transferListOfGroupOfCellFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    const int rank = getMPIRank();
    const int size = getMPISize();

    auto &otherRepe = other->getListOfGroupOfElementsbyElement();
    auto &otherLiel = other->getListOfGroupOfElements();

    const auto nbCells = connectionMesh->getNumberOfCells();
    auto &cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto &cellsOwner = connectionMesh->getCellsOwner();

    int nbCellsLoc = 0;
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            nbCellsLoc++;
    }

    VectorLong typeCellFE;
    typeCellFE.reserve( nbCellsLoc );
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank ) {
            auto cellId = ( *cellsLocNum )[i] - 1;
            auto numGrel = ( *otherRepe )[2 * cellId];
            if ( numGrel > 0 ) {
                auto &grel = ( *otherLiel )[numGrel];
                auto typeFE = grel[grel.size() - 1];
                typeCellFE.push_back( typeFE );
            } else {
                typeCellFE.push_back( 0 );
            }
        }
    }

    std::vector< VectorLong > typeFEGathered;
    AsterMPI::all_gather( typeCellFE, typeFEGathered );
    typeCellFE.clear();

    std::vector< VectorLong > listOfGrel;
    std::map< ASTERINTEGER, ASTERINTEGER > listOfGrelNume, listOfGrelNumeInv;

    _groupsOfCellsNumberByElement->allocate( 2 * nbCells );

    VectorLong nbCellsProc( size, 0 );
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        auto typeEF = typeFEGathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;

        if ( typeEF > 0 ) {
            if ( listOfGrelNume.count( typeEF ) == 0 ) {
                listOfGrelNume[typeEF] = listOfGrelNume.size();
                listOfGrelNumeInv[listOfGrelNumeInv.size()] = typeEF;
                listOfGrel.push_back( VectorLong() );
                listOfGrel[listOfGrelNume[typeEF]].reserve( nbCells );
            }

            listOfGrel[listOfGrelNume[typeEF]].push_back( i + 1 );

            ( *_groupsOfCellsNumberByElement )[2 * i] = listOfGrelNume[typeEF] + 1;
            ( *_groupsOfCellsNumberByElement )[2 * i + 1] =
                listOfGrel[listOfGrelNume[typeEF]].size();
        } else {
            ( *_groupsOfCellsNumberByElement )[2 * i] = 0;
            ( *_groupsOfCellsNumberByElement )[2 * i + 1] = 0;
        }
    }

    int totalSize = listOfGrelNume.size();
    int pos = 0;
    for ( auto &grel : listOfGrel ) {
        // il faut rajouter le nom de l'EF à la fin
        grel.push_back( listOfGrelNumeInv[pos++] );
        totalSize += grel.size();
    }

    _listOfGroupOfCells->allocateContiguousNumbered( listOfGrel.size(), totalSize );
    for ( auto &grel : listOfGrel ) {
        _listOfGroupOfCells->allocateObject( grel );
    }
};

void FiniteElementDescriptor::setFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    // Fill '.LGRF'
    _parameters->allocate( 3 );
    ( *_parameters )[0] = getMesh()->getName();
    _parameters->setInformationParameter( other->getParameters()->getInformationParameter() );

    // Fill '.NBNO'
    _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = 0;

    // Fill '.PRNM'
    transferDofDescriptorFrom( other );

    // Fill 'LIEL'
    transferListOfGroupOfCellFrom( other );
}

#endif /* ASTER_HAVE_MPI */
