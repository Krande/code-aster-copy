/**
 * @file ModelInterface.cxx
 * @brief Interface python de Model
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

#include "PythonBindings/ModelInterface.h"

#include "aster_pybind.h"

void exportModelToPython( py::module_ &mod ) {

    py::enum_< ModelSplitingMethod >( mod, "ModelSplitingMethod", R"(
Enumeration for model split method .
    )" )
        .value( "Centralized", Centralized )
        .value( "SubDomain", SubDomain )
        .value( "GroupOfCells", GroupOfCellsSplit )
        .export_values();

    py::enum_< GraphPartitioner >( mod, "GraphPartitioner", R"(
Enumeration for graph partitionner.
    )" )
        .value( "Scotch", ScotchPartitioner )
        .value( "Metis", MetisPartitioner )
        .export_values();

    py::class_< Model, Model::ModelPtr, DataStructure >( mod, "Model" )
#ifdef ASTER_HAVE_MPI
        .def( py::init( &initFactoryPtr< Model, ConnectionMeshPtr > ) )
        .def( py::init( &initFactoryPtr< Model, std::string, ConnectionMeshPtr > ) )
        .def( "getConnectionMesh", &Model::getConnectionMesh, R"(
            Return the ConnectionMesh

            Returns:
                ConnectionMesh: a pointer to the ConnectionMesh
            )" )
#endif /* ASTER_HAVE_MPI */
        .def( py::init( &initFactoryPtr< Model, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< Model, BaseMeshPtr, bool > ) )
        .def( py::init( &initFactoryPtr< Model, std::string, FiniteElementDescriptorPtr > ) )
        .def( py::init( &initFactoryPtr< Model, std::string, FiniteElementDescriptorPtr, bool > ) )
        .def( "addModelingOnMesh", &Model::addModelingOnMesh, R"(
Add modeling on all mesh

Arguments:
    physics (Physics): Physics
    modeling (Modelings): Modeling
    formulation (Formulation): Formulation (optional)
        )",
              py::arg( "physics" ), py::arg( "modeling" ),
              py::arg( "formulation" ) = NoFormulation )
        .def( "addModelingOnGroupOfCells", &Model::addModelingOnGroupOfCells, R"(
Add modeling on all mesh

Arguments:
    physics (Physics): Physics
    modeling (Modelings): Modeling
    grpma (str): Name of element group
    formulation (Formulation): Formulation (optional)
        )",
              py::arg( "physics" ), py::arg( "modeling" ), py::arg( "grpma" ),
              py::arg( "formulation" ) = NoFormulation )
        .def( "build", &Model::build )
        .def( "existsThm", &Model::existsThm )
        .def( "existsHHO", &Model::existsHHO )
        .def( "isAxis", &Model::isAxis, R"(
            To know if the model is Axisymmetric

            Returns:
                bool: *True* if the model is axisymmetric, else *False*
            )" )
        .def( "existsRdM", &Model::existsRdM, R"(
            To know if the model has RdM elements

            Returns:
                bool: *True* if the model contains beam, shell or discret elements, else *False*
            )" )
        .def( "existsPartition", &Model::existsPartition )
        .def( "getModelisationName", &Model::getModelisationName, R"(
            Get modelisation name used in model

            Returns:
                str: modelisation name if single modelisation, else '#PLUSIEURS'
            )" )
        .def( "getPartitionMethod", &Model::getPartitionMethod, R"(
            Get partition method

            Returns:
                str: partition method
            )" )
        .def( "isXfem", &Model::isXfem )
        .def( "getXfemContact", &Model::getXfemContact )
        .def( "existsMultiFiberBeam", &Model::existsMultiFiberBeam )
        .def( "getSaneModel", &Model::getSaneModel )
        .def( "getMesh", &Model::getMesh, R"(
            Return the mesh

            Returns:
                Mesh: a pointer to the mesh
            )" )
        .def( "isMechanical", &Model::isMechanical, R"(
            To know if the model is mechanical or not

            Returns:
                bool: True - if the model is mechanical
            )" )
        .def( "isThermal", &Model::isThermal, R"(
            To know if the model is thermal or not

            Returns:
                bool: True - if the model is thermal
            )" )
        .def( "isAcoustic", &Model::isAcoustic, R"(
            To know if the model is acoustic or not

            Returns:
                bool: True - if the model is acoustic
            )" )
        .def( "getPhysics", &Model::getPhysics, R"(
            To know the physics supported by the model

            Returns:
                str: Mechanics or Thermal or Acoustic
            )" )
        .def( "getGeometricDimension", &Model::getGeometricDimension, R"(
            To know the geometric dimension supported by the model

            Returns:
                int: geometric dimension
            )" )
        .def( "getSplittingMethod", &Model::getSplittingMethod )
        .def( "getGraphPartitioner", &Model::getGraphPartitioner )
        .def( "setSaneModel", &Model::setSaneModel )
        .def( "xfemPreconditioningEnable", &Model::xfemPreconditioningEnable )
        .def( "banBalancing", &Model::banBalancing, R"(
            Prohibit model balancing
            )" )
        .def( "setSplittingMethod", py::overload_cast< ModelSplitingMethod, GraphPartitioner >(
                                        &Model::setSplittingMethod ) )
        .def( "setSplittingMethod",
              py::overload_cast< ModelSplitingMethod >( &Model::setSplittingMethod ) )
        .def( "getFiniteElementDescriptor", &Model::getFiniteElementDescriptor )
#ifdef ASTER_HAVE_MPI
        .def( "setFrom", &Model::setFrom, R"(
            Set a model defined on a ConnectionMesh from an other model

            Arguments:
                model (Model): Table identifier.
            )",
              py::arg( "model" ) )
#endif
        .def( "getTable", &ListOfTables::getTable, R"(
            Extract a Table from the datastructure.

            Arguments:
                identifier (str): Table identifier.

            Returns:
                Table: Table stored with the given identifier.
            )",
              py::arg( "identifier" ) );
};
