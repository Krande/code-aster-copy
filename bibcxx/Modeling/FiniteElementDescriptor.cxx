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
#include "ParallelUtilities/MPIInfos.h"
#include "ParallelUtilities/AsterMPI.h"

FiniteElementDescriptorClass::FiniteElementDescriptorClass( const std::string &name,
                                                            const BaseMeshPtr mesh,
                                                            const JeveuxMemory memType )
    : DataStructure( name, 19, "LIGREL", memType ),
      _numberOfDelayedNumberedConstraintNodes( getName() + ".NBNO" ),
      _parameters( getName() + ".LGRF" ), _dofDescriptor( getName() + ".PRNM" ),
      _listOfGroupOfCells( getName() + ".LIEL" ),
      _groupsOfCellsNumberByElement( getName() + ".REPE" ),
      _delayedNumberedConstraintElementsDescriptor( getName() + ".NEMA" ),
      _dofOfDelayedNumberedConstraintNodes( getName() + ".PRNS" ),
      _delayedNodesNumbering( getName() + ".LGNS" ),
      _superElementsDescriptor( getName() + ".SSSA" ),
      _nameOfNeighborhoodStructure( getName() + ".NVGE" ), _mesh( mesh ),
      _explorer(
          ConnectivityDelayedElementsExplorer( _delayedNumberedConstraintElementsDescriptor ) ),
      _explorer2( ConnectivityDelayedElementsExplorer( _listOfGroupOfCells ) ){};

int FiniteElementDescriptorClass::getPhysics( void ) const
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
void FiniteElementDescriptorClass::transferDofDescriptorFrom( FiniteElementDescriptorPtr &other ) {
    if ( !getMesh()->isPartial() )
        throw std::runtime_error(
            "the mesh associated to finiteElementDescriptorClass is not a partial mesh" );
    const ConnectionMeshPtr connectionMesh =
        boost::static_pointer_cast< ConnectionMeshClass >( getMesh() );
    if ( connectionMesh->getParallelMesh() != other->getMesh() )
        throw std::runtime_error(
            "parallel mesh associated to partial mesh of FiniteElementDescriptorClass \n"
            "does not correspond to other FiniteElementDescriptorClass mesh" );
    getPhysicalNodesComponentDescriptor();

    const int rank = getMPIRank();
    int nbNodes = connectionMesh->getNumberOfNodes();
    int nec = _dofDescriptor->size() / nbNodes;

    const JeveuxVectorLong &localNumbering = connectionMesh->getLocalNumbering();
    const JeveuxVectorLong &owner = connectionMesh->getOwner();
    const JeveuxVectorLong &otherDofDescriptor = other->getPhysicalNodesComponentDescriptor();

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

    VectorLong bufferGathered;
    AsterMPI::all_gather( buffer, _dofDescriptor );
    buffer.clear();

};
#endif /* ASTER_HAVE_MPI */
