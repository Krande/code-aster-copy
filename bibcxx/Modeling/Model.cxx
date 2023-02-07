/**
 * @file Model.cxx
 * @brief Implementation de Model
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

#include "Modeling/Model.h"

#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "astercxx.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>
#include <typeinfo>

const char *const ModelSplitingMethodNames[nbModelSplitingMethod] = {"CENTRALISE", "SOUS_DOMAINE",
                                                                     "GROUP_ELEM"};
const char *const GraphPartitionerNames[nbGraphPartitioner] = {"SCOTCH", "METIS"};

Model::Model( const std::string name, const bool is_xfem )
    : DataStructure( name, 8, "MODELE" ),
      ListOfTables( name ),
      _typeOfCells( JeveuxVectorLong( getName() + ".MAILLE    " ) ),
      _typeOfNodes( JeveuxVectorLong( getName() + ".NOEUD     " ) ),
      _baseMesh( nullptr ),
      _splitMethod( SubDomain ),
      _graphPartitioner( MetisPartitioner ),
      _ligrel( nullptr ),
      _xfemModel( nullptr ) {
#ifdef ASTER_HAVE_MPI
    _connectionMesh = nullptr;
#endif

    if ( is_xfem )
        _xfemModel = std::make_shared< XfemModel >( getName() );
};

Model::Model( const std::string name, const FiniteElementDescriptorPtr fed, const bool is_xfem )
    : Model( name, is_xfem ) {
    _baseMesh = fed->getMesh();
    _ligrel = fed;

    AS_ASSERT( !_baseMesh->isEmpty() );
};

Model::Model( const BaseMeshPtr mesh, const bool is_xfem )
    : Model( DataStructureNaming::getNewName(), is_xfem ) {
    _baseMesh = mesh;
    _ligrel = std::make_shared< FiniteElementDescriptor >( getName() + ".MODELE", _baseMesh );

    AS_ASSERT( !_baseMesh->isEmpty() );
};

#ifdef ASTER_HAVE_MPI
Model::Model( const std::string name, const ConnectionMeshPtr mesh ) : Model( name ) {
    _baseMesh = mesh;
    _connectionMesh = mesh;
    _ligrel = std::make_shared< FiniteElementDescriptor >( getName() + ".MODELE", _baseMesh );

    AS_ASSERT( !_baseMesh->isEmpty() );
    AS_ASSERT( !_connectionMesh->isEmpty() );
};

Model::Model( const ConnectionMeshPtr mesh ) : Model( DataStructureNaming::getNewName(), mesh ){};
#endif /* ASTER_HAVE_MPI */

bool Model::existsSuperElement() { return ( numberOfSuperElement() > 0 ); }

FiniteElementDescriptorPtr Model::getFiniteElementDescriptor() const { return _ligrel; };

SyntaxMapContainer Model::buildModelingsSyntaxMapContainer() const {
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

bool Model::buildWithSyntax( SyntaxMapContainer &dict ) {
    CommandSyntax cmdSt( "AFFE_MODELE" );
    cmdSt.setResult( ResultNaming::getCurrentName(), "MODELE" );
    cmdSt.define( dict );

    // Maintenant que le fichier de commande est pret, on appelle OP0018
    ASTERINTEGER op = 18;
    CALL_EXECOP( &op );

    // Attention, la connection des objets a leur image JEVEUX n'est pas necessaire
    _typeOfCells->updateValuePointer();

    return true;
};

bool Model::build() {
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

    return buildWithSyntax( dict ) && update_tables() && _ligrel->build();
};

bool Model::existsThm() {
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

bool Model::existsMultiFiberBeam() {
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

bool Model::xfemPreconditioningEnable() {
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

bool Model::existsXfem() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_XFEM" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool Model::existsHHO() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_HHO" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool Model::existsAxis() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_AXIS" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool Model::exists3DShell() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_COQ3D" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool Model::existsSTRX() {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_STRX" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = trim( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

ASTERINTEGER Model::numberOfSuperElement() { return _ligrel->numberOfSuperElement(); };

bool Model::existsFiniteElement() { return _ligrel->existsFiniteElement(); };

void Model::addModelingOnMesh( Physics phys, Modelings mod ) {
    _modelisations.push_back( listOfModsAndGrpsValue( ElementaryModeling( phys, mod ),
                                                      MeshEntityPtr( new AllMeshEntities() ) ) );
};

void Model::addModelingOnGroupOfCells( Physics phys, Modelings mod, std::string nameOfGroup ) {
    if ( !_baseMesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_baseMesh->hasGroupOfCells( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _modelisations.push_back( listOfModsAndGrpsValue(
        ElementaryModeling( phys, mod ), MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) ) );
};

void Model::addModelingOnGroupOfNodes( Physics phys, Modelings mod, std::string nameOfGroup ) {
    if ( !_baseMesh )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !_baseMesh->hasGroupOfNodes( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _modelisations.push_back( listOfModsAndGrpsValue(
        ElementaryModeling( phys, mod ), MeshEntityPtr( new GroupOfNodes( nameOfGroup ) ) ) );
};

void Model::setXfemModel() {
    const std::string modelName = getName();
    _xfemModel = std::make_shared< XfemModel >( modelName );
};

int Model::getGeometricDimension( void ) const {
    ASTERINTEGER nb_dim_;
    ASTERINTEGER ier;
    std::string repk;
    CALL_DISMOI( "DIM_GEOM", this->getName().c_str(), "MODELE", &nb_dim_, repk.c_str(), "F", &ier );
    return nb_dim_;
}

GraphPartitioner Model::getGraphPartitioner() const { return _graphPartitioner; };

#ifdef ASTER_HAVE_MPI
ConnectionMeshPtr Model::getConnectionMesh() const {
    if ( ( !_connectionMesh ) || _connectionMesh->isEmpty() )
        throw std::runtime_error( "ConnectionMesh of model is empty" );
    return _connectionMesh;
};
#endif /* ASTER_HAVE_MPI */

ModelPtr Model::getSaneModel() const { return _saneModel; };

XfemModelPtr Model::getXfemModel() const { return _xfemModel; };

ModelSplitingMethod Model::getSplittingMethod() const { return _splitMethod; };

BaseMeshPtr Model::getMesh() const {
    if ( ( !_baseMesh ) || _baseMesh->isEmpty() )
        throw std::runtime_error( "Mesh of model is empty" );
    return _baseMesh;
};

JeveuxVectorLong Model::getTypeOfCells() const { return _typeOfCells; }

bool Model::isEmpty() const { return !_typeOfCells->exists(); };

void Model::setSaneModel( ModelPtr saneModel ) { _saneModel = saneModel; };

void Model::setSplittingMethod( ModelSplitingMethod split, GraphPartitioner partitioner ) {
    setSplittingMethod( split );
    _graphPartitioner = partitioner;
};

void Model::setSplittingMethod( ModelSplitingMethod split ) {
#ifdef ASTER_HAVE_MPI
    if ( _connectionMesh && !_connectionMesh->isEmpty() && split != Centralized )
        throw std::runtime_error( "For Parallel mesh, Centralized splitting is mandatory" );
#endif /* ASTER_HAVE_MPI */

    _splitMethod = split;
};

bool Model::isMechanical( void ) const { return this->getPhysics() == Physics::Mechanics; };

bool Model::isThermal( void ) const { return this->getPhysics() == Physics::Thermal; };

bool Model::isAcoustic( void ) const { return this->getPhysics() == Physics::Acoustic; };

int Model::getPhysics( void ) const { return _ligrel->getPhysics(); }

bool Model::isXfem() const { return _xfemModel != nullptr; };

#ifdef ASTER_HAVE_MPI
bool Model::setFrom( const ModelPtr model ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == model->getMesh() );

    // tranfer LIGREL
    auto Fed = model->getFiniteElementDescriptor();
    this->getFiniteElementDescriptor()->setFrom( Fed );

    // tranfer .MAILLE
    const auto nbCells = connectionMesh->getNumberOfCells();
    _typeOfCells->allocate( nbCells );
    auto &cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto &cellsOwner = connectionMesh->getCellsOwner();

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
            buffer.push_back( ( *typeCellsOther )[( *cellsLocNum )[i] - 1] );
    }

    std::vector< VectorLong > gathered;
    AsterMPI::all_gather( buffer, gathered );

    VectorLong nbCellsProc( size, 0 );
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        ( *_typeOfCells )[i] = gathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;
    }

    return true;
};
#endif
