/**
 * @file FiniteElementDescriptor.cxx
 * @brief Implementation de FiniteElementDescriptor
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <algorithm>

#include "Modeling/FiniteElementDescriptor.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/ConnectionMesh.h"
#include "Modeling/PhysicalQuantityManager.h"
#include "Modeling/PhysicsAndModelings.h"
#include "ParallelUtilities/AsterMPI.h"

FiniteElementDescriptor::FiniteElementDescriptor( const std::string &name,
                                                            const BaseMeshPtr mesh,
                                                            const JeveuxMemory memType )
    : DataStructure( name, 19, "LIGREL", memType ),
      _numberOfDelayedNumberedConstraintNodes( getName() + ".NBNO" ),
      _parameters( getName() + ".LGRF" ), _dofDescriptor( getName() + ".PRNM" ),
      _listOfGroupOfCells( getName() + ".LIEL" ),
      _groupsOfCellsNumberByElement( getName() + ".REPE" ),
      _delayedNumberedConstraintElementsDescriptor( getName() + ".NEMA" ),
      _dofOfDelayedNumberedConstraintNodes( getName() + ".PRNS" ),
      _virtualNodesNumbering( getName() + ".LGNS" ),
      _superElementsDescriptor( getName() + ".SSSA" ),
      _nameOfNeighborhoodStructure( getName() + ".NVGE" ), _mesh( mesh ),
      _explorer(
          ConnectivityVirtualCellsExplorer( _delayedNumberedConstraintElementsDescriptor ) ),
      _explorer2( ConnectivityVirtualCellsExplorer( _listOfGroupOfCells ) ){};

int FiniteElementDescriptor::getPhysics( void ) const
{
    const std::string docu = trim(_parameters->getInformationParameter());

    if( docu == "MECA" )
        return Physics::Mechanics;
    else if( docu == "THER" )
        return Physics::Thermal;
    else if( docu == "ACOU" )
        return Physics::Acoustic;
    else
        throw std::runtime_error("Unknown physics");

    return -1;
};

#ifdef ASTER_HAVE_MPI
void FiniteElementDescriptor::transferDofDescriptorFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        boost::static_pointer_cast< ConnectionMesh >( getMesh() );

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

    std::vector<VectorLong> gathered;
    AsterMPI::all_gather( buffer, gathered );
    buffer.clear();

    _dofDescriptor->allocate(nbNodes * nec);
    VectorLong nbNodesProc( size, 0);
    for ( int i = 0; i < nbNodes; ++i ) {
        auto rowner = ( *owner )[i];
        for( int j = 0; j < nec; ++j )
            (*_dofDescriptor)[i*nec + j] = gathered[rowner][nbNodesProc[rowner]*nec + j];
        nbNodesProc[rowner]++;
    }

};

void FiniteElementDescriptor::transferListOfGroupOfCellFrom( FiniteElementDescriptorPtr& other)
{
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        boost::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    const int rank = getMPIRank();
    const int size = getMPISize();

    auto& otherRepe = other->getListOfGroupOfCellsbyCell();
    auto& otherLiel = other->getListOfGroupOfCells();

    const auto nbCells = connectionMesh->getNumberOfCells();
    auto& cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto& cellsOwner = connectionMesh->getCellsOwner();

    int nbCellsLoc = 0;
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            nbCellsLoc++;
    }

    VectorLong typeCellFE;
    typeCellFE.reserve( nbCellsLoc );
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
        {
            auto cellId = (*cellsLocNum)[i] - 1;
            auto numGrel = (*otherRepe)[ 2*cellId];
            if( numGrel > 0)
            {
                auto& grel = otherLiel->getObject( numGrel );
                auto typeFE = grel[grel.size()-1];
                typeCellFE.push_back( typeFE );
            }
            else
            {
                typeCellFE.push_back( 0 );
            }
        }
    }

    std::vector<VectorLong> typeFEGathered;
    AsterMPI::all_gather(typeCellFE, typeFEGathered);
    typeCellFE.clear();

    std::vector< VectorLong > listOfGrel;
    std::map< ASTERINTEGER, ASTERINTEGER> listOfGrelNume, listOfGrelNumeInv;

    _groupsOfCellsNumberByElement->allocate( 2*nbCells );

    VectorLong nbCellsProc( size, 0);
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        auto typeEF = typeFEGathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;

        if( typeEF > 0)
        {
            if(listOfGrelNume.count(typeEF) == 0 )
            {
                listOfGrelNume[typeEF] = listOfGrelNume.size();
                listOfGrelNumeInv[listOfGrelNumeInv.size()] = typeEF;
                listOfGrel.push_back( VectorLong() );
                listOfGrel[listOfGrelNume[typeEF]].reserve(nbCells);
            }

            listOfGrel[listOfGrelNume[typeEF]].push_back(i+1);

            (*_groupsOfCellsNumberByElement)[2*i] = listOfGrelNume[typeEF] + 1;
            (*_groupsOfCellsNumberByElement)[2*i + 1] = listOfGrel[listOfGrelNume[typeEF]].size();
        }
        else
        {
            (*_groupsOfCellsNumberByElement)[2*i] = 0;
            (*_groupsOfCellsNumberByElement)[2*i + 1] = 0;
        }
    }

    int totalSize = listOfGrelNume.size();
    int pos = 0;
    for( auto& grel : listOfGrel)
    {
        // il faut rajouter le nom de l'EF à la fin
        grel.push_back(listOfGrelNumeInv[pos++]);
        totalSize += grel.size();
    }

    _listOfGroupOfCells->allocateContiguous(Permanent, listOfGrel.size(), totalSize, Numbered);
    int posInCollection = 1;
    for( auto& grel : listOfGrel)
    {
        _listOfGroupOfCells->allocateObject( grel.size() );
        _listOfGroupOfCells->getObject( posInCollection ).setValues( grel );
        ++posInCollection;
    }
};


void FiniteElementDescriptor::setFrom( FiniteElementDescriptorPtr &other )
{
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        boost::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    // Fill '.LGRF'
    _parameters->allocate(3);
    (*_parameters)[0] = getMesh()->getName();
    _parameters->setInformationParameter(other->getParameters()->getInformationParameter());

    // Fill '.NBNO'
    _numberOfDelayedNumberedConstraintNodes->allocate(1);
    (*_numberOfDelayedNumberedConstraintNodes)[0] = 0;

    // Fill '.PRNM'
    transferDofDescriptorFrom(other);

    // Fill 'LIEL'
    transferListOfGroupOfCellFrom(other);
}
#endif /* ASTER_HAVE_MPI */
