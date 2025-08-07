/**
 * @file Model.cxx
 * @brief Implementation de Model
 * @author Nicolas Sellenet
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

#include "astercxx.h"

#include "Modeling/Model.h"

#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>
#include <typeinfo>

Model::Model( const std::string name, const bool is_xfem )
    : DataStructure( name, 8, "MODELE" ),
      ListOfTables( name ),
      _typeOfNodes( getName() + ".NOEUD     " ),
      _splitMethod( SubDomain ),
      _graphPartitioner( MetisPartitioner ),
      _ligrel( nullptr ),
      _hhoModel( std::make_shared< HHOModel >( getName() ) ),
      _xfemModel( nullptr ) {
#ifdef ASTER_HAVE_MPI
    _connectionMesh = nullptr;
#endif

    if ( is_xfem )
        _xfemModel = std::make_shared< XfemModel >( getName() );
};

Model::Model( const std::string name, const FiniteElementDescriptorPtr fed, const bool is_xfem )
    : Model( name, is_xfem ) {
    _ligrel = fed;

    AS_ASSERT( !_ligrel->getMesh()->isEmpty() );
};

Model::Model( const BaseMeshPtr mesh, const bool is_xfem )
    : Model( DataStructureNaming::getNewName(), is_xfem ) {
    _ligrel = std::make_shared< FiniteElementDescriptor >( getName() + ".MODELE", mesh );

    AS_ASSERT( !mesh->isEmpty() );
};

#ifdef ASTER_HAVE_MPI
Model::Model( const std::string name, const ConnectionMeshPtr mesh ) : Model( name ) {
    _connectionMesh = mesh;
    _ligrel = std::make_shared< FiniteElementDescriptor >( getName() + ".MODELE", mesh );

    AS_ASSERT( !mesh->isEmpty() );
    AS_ASSERT( !_connectionMesh->isEmpty() );
};

Model::Model( const ConnectionMeshPtr mesh ) : Model( DataStructureNaming::getNewName(), mesh ) {};
#endif /* ASTER_HAVE_MPI */

Model::Model( const BaseMeshPtr mesh, const ModelPtr model ) : Model( mesh ) {
    _modelisations = model->_modelisations;
};

bool Model::existsSuperElement() { return ( numberOfSuperElement() > 0 ); }

FiniteElementDescriptorPtr Model::getFiniteElementDescriptor() const { return _ligrel; };

SyntaxMapContainer Model::buildModelingsSyntaxMapContainer() const {
    SyntaxMapContainer dict;

    dict.container["VERI_JACOBIEN"] = "OUI";
    if ( !this->getMesh() )
        throw std::runtime_error( "Mesh is undefined" );
    dict.container["MAILLAGE"] = this->getMesh()->getName();

    ListSyntaxMapContainer listeAFFE;
    for ( listOfModsAndGrpsCIter curIter = _modelisations.begin(); curIter != _modelisations.end();
          ++curIter ) {
        SyntaxMapContainer dict2;
        dict2.container["PHENOMENE"] = curIter->first.getPhysic();
        dict2.container["MODELISATION"] = curIter->first.getModeling();

        VectorString grpMa;
        const auto &curList = curIter->second;
        for ( const auto &curIter : curList ) {
            if ( ( *( curIter ) ).getType() == AllMeshEntitiesType ) {
                dict2.container["TOUT"] = "OUI";
            } else {
                if ( ( *( curIter ) ).getType() == GroupOfCellsType ) {
                    grpMa.push_back( ( curIter )->getName() );
                }
            }
        }
        if ( grpMa.size() != 0 )
            dict2.container["GROUP_MA"] = grpMa;
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

    return true;
};

bool Model::build() {
    SyntaxMapContainer dict = buildModelingsSyntaxMapContainer();
    if ( this->getMesh()->isParallel() ) {
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

bool Model::existsThm() const { return dismoi( "EXI_THM", false ) == "OUI"; };

bool Model::existsMultiFiberBeam() const { return dismoi( "EXI_STR2", false ) == "OUI"; };

bool Model::xfemPreconditioningEnable() const { return dismoi( "PRE_COND_XFEM", false ) == "OUI"; };

bool Model::existsXfem() const { return dismoi( "EXI_XFEM", false ) == "OUI"; };

bool Model::existsHHO() const { return dismoi( "EXI_HHO", false ) == "OUI"; };

bool Model::existsAxis() const { return dismoi( "EXI_AXIS", false ) == "OUI"; };

bool Model::isAxis() const { return dismoi( "AXIS" ) == "OUI"; };

bool Model::exists3DShell() const { return dismoi( "EXI_COQ3D" ) == "OUI"; };

bool Model::existsSTRX() const { return dismoi( "EXI_STRX" ) == "OUI"; };

bool Model::existsRdM() const { return dismoi( "EXI_RDM" ) == "OUI"; };

bool Model::existsPartition() const { return dismoi( "PARTITION" ) != ""; }

const std::string Model::getModelisationName() const { return dismoi( "MODELISATION" ); };

const std::string Model::getPartitionMethod() const { return _ligrel->getPartitionMethod(); };

const std::string Model::dismoi( const std::string &question, bool stop ) const {
    const std::string typeco( "MODELE" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    std::string arret;
    if ( stop )
        arret = "F";
    else
        arret = "C";
    CALLO_DISMOI( question, getName(), typeco, &repi, repk, arret, &ier );
    return strip( repk.toString() );
};

ASTERINTEGER Model::numberOfSuperElement() { return _ligrel->numberOfSuperElement(); };

bool Model::existsFiniteElement() { return _ligrel->existsFiniteElement(); };

void Model::addModelingOnMesh( Physics phys, Modelings mod, Formulation form ) {
    _modelisations.push_back( listOfModsAndGrpsValue(
        ElementaryModeling( phys, mod, form ), { MeshEntityPtr( new AllMeshEntities() ) } ) );
};

void Model::addModelingOnGroupOfCells( Physics phys, Modelings mod, std::string nameOfGroup,
                                       Formulation form ) {
    if ( !this->getMesh() )
        throw std::runtime_error( "Mesh is not defined" );
    if ( !this->getMesh()->hasGroupOfCells( nameOfGroup ) )
        throw std::runtime_error( nameOfGroup + " not in mesh" );

    _modelisations.push_back(
        listOfModsAndGrpsValue( ElementaryModeling( phys, mod, form ),
                                { MeshEntityPtr( new GroupOfCells( nameOfGroup ) ) } ) );
};

void Model::setXfemModel() {
    const std::string modelName = getName();
    _xfemModel = std::make_shared< XfemModel >( modelName );
};

/**
 * @brief Get the HHO model
 */
HHOModelPtr Model::getHHOModel() const { return _hhoModel; };

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
    if ( ( !_ligrel ) || ( !_ligrel->getMesh() ) || _ligrel->getMesh()->isEmpty() )
        throw std::runtime_error( "Mesh of model is empty" );
    return _ligrel->getMesh();
};

JeveuxVectorLong Model::getFiniteElementType() const { return _ligrel->getFiniteElementType(); }

bool Model::isEmpty() const { return !_ligrel->getFiniteElementType().exists(); };

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

ASTERINTEGER Model::getXfemContact() const { return isXfem() ? _xfemModel->getContact() : -1; };

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

    return true;
};
#endif
