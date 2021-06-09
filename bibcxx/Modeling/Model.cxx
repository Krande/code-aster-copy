/**
 * @file Model.cxx
 * @brief Implementation de ModelClass
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <stdexcept>
#include <typeinfo>

#include "astercxx.h"

#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Modeling/Model.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

const char *const ModelSplitingMethodNames[nbModelSplitingMethod] = {"CENTRALISE", "SOUS_DOMAINE",
                                                                     "GROUP_ELEM"};
const char *const GraphPartitionerNames[nbGraphPartitioner] = {"SCOTCH", "METIS"};

SyntaxMapContainer ModelClass::buildModelingsSyntaxMapContainer() const {
    SyntaxMapContainer dict;

    dict.container["VERI_JACOBIEN"] = "OUI";
    if ( !_baseMesh )
        throw std::runtime_error( "Mesh is undefined" );
    dict.container["MAILLAGE"] = _baseMesh->getName();

    ListSyntaxMapContainer listeAFFE;
    for ( listOfModsAndGrpsCIter curIter = _modelisations.begin(); curIter != _modelisations.end();
          ++curIter ) {
        SyntaxMapContainer dict2;
        dict2.container["PHENOMENE"] = curIter->first.getPhysic();
        dict2.container["MODELISATION"] = curIter->first.getModeling();

        if ( ( *( curIter->second ) ).getType() == AllMeshEntitiesType ) {
            dict2.container["TOUT"] = "OUI";
        } else {
            if ( ( *( curIter->second ) ).getType() == GroupOfNodesType )
                dict2.container["GROUP_NO"] = ( curIter->second )->getName();
            else if ( ( *( curIter->second ) ).getType() == GroupOfCellsType )
                dict2.container["GROUP_MA"] = ( curIter->second )->getName();
        }
        listeAFFE.push_back( dict2 );
    }
    dict.container["AFFE"] = listeAFFE;
    return dict;
};

bool ModelClass::buildWithSyntax( SyntaxMapContainer &dict ) {
    CommandSyntax cmdSt( "AFFE_MODELE" );
    cmdSt.setResult( ResultNaming::getCurrentName(), "MODELE" );
    cmdSt.define( dict );

    // Maintenant que le fichier de commande est pret, on appelle OP0018
    try {
        ASTERINTEGER op = 18;
        CALL_EXECOP( &op );
    } catch ( ... ) {
        throw;
    }
    // Attention, la connection des objets a leur image JEVEUX n'est pas necessaire
    _typeOfCells->updateValuePointer();

    return true;
};

bool ModelClass::build() {
    SyntaxMapContainer dict = buildModelingsSyntaxMapContainer();
    if ( _baseMesh->isParallel() ) {
        ListSyntaxMapContainer listeDISTRIBUTION;
        SyntaxMapContainer dict2;
        dict2.container["METHODE"] = ModelSplitingMethodNames[(int)Centralized];
        listeDISTRIBUTION.push_back( dict2 );
        dict.container["DISTRIBUTION"] = listeDISTRIBUTION;
    } else {
        ListSyntaxMapContainer listeDISTRIBUTION;
        SyntaxMapContainer dict2;
        dict2.container["METHODE"] = ModelSplitingMethodNames[(int)getSplittingMethod()];
        dict2.container["PARTITIONNEUR"] = GraphPartitionerNames[(int)getGraphPartitioner()];
        listeDISTRIBUTION.push_back( dict2 );
        dict.container["DISTRIBUTION"] = listeDISTRIBUTION;
    }

    return buildWithSyntax( dict ) && update_tables();
};

bool ModelClass::existsThm() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_THM" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool ModelClass::existsMultiFiberBeam() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_STR2" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool ModelClass::xfemPreconditioningEnable() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "PRE_COND_XFEM" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};


void ModelClass::addModelingOnGroupOfCells( Physics phys, Modelings mod, std::string nameOfGroup ) {
    if ( !_baseMesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_baseMesh->hasGroupOfCells( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _modelisations.push_back( listOfModsAndGrpsValue(
        ElementaryModeling( phys, mod ), MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) ) );
};

void ModelClass::addModelingOnGroupOfNodes( Physics phys, Modelings mod, std::string nameOfGroup ) {
    if ( !_baseMesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_baseMesh->hasGroupOfNodes( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _modelisations.push_back( listOfModsAndGrpsValue(
        ElementaryModeling( phys, mod ), MeshEntityPtr( new GroupOfNodes( nameOfGroup ) ) ) );
    };


#ifdef ASTER_HAVE_MPI
bool ModelClass::transferFrom( const ModelPtr model)
{
    // "the mesh associated to finiteElementDescriptorClass is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        boost::static_pointer_cast< ConnectionMeshClass >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptorClass \n"
    //        "does not correspond to other FiniteElementDescriptorClass mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == model->getMesh() );

    // tranfer LIGREL
    auto Fed = model->getFiniteElementDescriptor();
    this->getFiniteElementDescriptor()->transferFrom(Fed);

    // tranfer .MAILLE
    const auto nbCells = connectionMesh->getNumberOfCells();
    _typeOfCells->allocate(nbCells);
    auto& cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto& cellsOwner = connectionMesh->getCellsOwner();

    const auto rank = getMPIRank();
    const auto size = getMPISize();

    int nbCellsLoc = 0;
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            nbCellsLoc++;
    }

    auto typeCellsOther = model->getTypeOfCells();
    typeCellsOther->updateValuePointer();

    VectorLong buffer;
    buffer.reserve( nbCellsLoc );
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            buffer.push_back( (*typeCellsOther)[ (*cellsLocNum)[i] - 1 ] );
    }

    std::vector<VectorLong> gathered;
    AsterMPI::all_gather(buffer, gathered);

    VectorLong nbCellsProc( size, 0);
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        (*_typeOfCells)[i] = gathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;
    }

    return true;
};

#endif
