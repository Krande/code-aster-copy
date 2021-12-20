/**
 * @file ModelInterface.cxx
 * @brief Interface python de Model
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

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/ModelInterface.h"

void exportModelToPython() {

    py::enum_< ModelSplitingMethod >( "ModelSplitingMethod" )
        .value( "Centralized", Centralized )
        .value( "SubDomain", SubDomain )
        .value( "GroupOfCells", GroupOfCellsSplit );

    py::enum_< GraphPartitioner >( "GraphPartitioner" ).value( "Scotch", ScotchPartitioner ).value(
        "Metis", MetisPartitioner );


    void ( Model::*split1 )( ModelSplitingMethod ) = &Model::setSplittingMethod;

    void ( Model::*split2 )( ModelSplitingMethod, GraphPartitioner ) =
        &Model::setSplittingMethod;

    py::class_< Model, Model::ModelPtr, py::bases< DataStructure > >( "Model",
                                                                                 py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< Model, BaseMeshPtr>))
        .def( "__init__", py::make_constructor(&initFactoryPtr< Model, std::string,
                                                                BaseMeshPtr>))
#ifdef ASTER_HAVE_MPI
        .def( "__init__", py::make_constructor(&initFactoryPtr< Model, ConnectionMeshPtr>))
        .def( "__init__", py::make_constructor(&initFactoryPtr< Model, std::string,
                                                                            ConnectionMeshPtr>))
#endif /* ASTER_HAVE_MPI */
        .def( "addModelingOnMesh", &Model::addModelingOnMesh )
        .def( "addModelingOnGroupOfCells", &Model::addModelingOnGroupOfCells )
        .def( "addModelingOnGroupOfNodes", &Model::addModelingOnGroupOfNodes )
        .def( "build", &Model::build )
        .def( "existsThm", &Model::existsThm )
        .def( "existsMultiFiberBeam", &Model::existsMultiFiberBeam )
        .def( "getSaneModel", &Model::getSaneModel )
        .def( "getMesh", &Model::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )",
              ( py::arg( "self" ) )  )
        .def( "isMechanical", &Model::isMechanical, R"(
To know if the model is mechanical or not

Returns:
    Bool: True - if the model is mechanical
        )",
              ( py::arg( "self" ) )  )
        .def( "isThermal", &Model::isThermal, R"(
To know if the model is thermal or not

Returns:
    Bool: True - if the model is thermal
        )",
              ( py::arg( "self" ) )  )
        .def( "isAcoustic", &Model::isAcoustic, R"(
To know if the model is acoustic or not

Returns:
    Bool: True - if the model is acoustic
        )",
              ( py::arg( "self" ) )  )
        .def( "getPhysics", &Model::getPhysics, R"(
To know the physics supported by the model

Returns:
    str: Mechanics or Thermal or Acoustic
        )",
              ( py::arg( "self" ) )  )
        .def( "getSplittingMethod", &Model::getSplittingMethod )
        .def( "getGraphPartitioner", &Model::getGraphPartitioner )
        .def( "setSaneModel", &Model::setSaneModel )
        .def( "xfemPreconditioningEnable", &Model::xfemPreconditioningEnable )
        .def( "setSplittingMethod", split1 )
        .def( "setSplittingMethod", split2 )
        .def( "getFiniteElementDescriptor", &Model::getFiniteElementDescriptor )
#ifdef ASTER_HAVE_MPI
        .def( "setFrom", &Model::setFrom, R"(
Set a model defined on a ConnectionMesh from an other model

Arguments:
    model (Model): Table identifier.

        )",
              ( py::arg( "self" ), py::arg( "model" ) ) )
#endif
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) )
        ;
};
