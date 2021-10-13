/**
 * @file PhysicalProblemInterface.cxx
 * @brief Interface python de PhysicalProblem
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
#include "PythonBindings/factory.h"
#include "PythonBindings/PhysicalProblemInterface.h"
#include "PythonBindings/LoadUtilities.h"

void exportPhysicalProblemToPython() {

    py::class_< PhysicalProblem, PhysicalProblemPtr > c1( "PhysicalProblem",
                                                                    py::no_init );
    // fake initFactoryPtr: not a DataStructure
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr >));
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr,
                                ElementaryCharacteristicsPtr >));
    c1.def( "getModel", &PhysicalProblem::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getMesh", &PhysicalProblem::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getMaterialField", &PhysicalProblem::getMaterialField, R"(
Return the material field

Returns:
    MaterialFieldPtr: a pointer to the material field
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getCodedMaterial", &PhysicalProblem::getCodedMaterial, R"(
Return the coded material

Returns:
    CodedMaterialPtr: a pointer to the coded material
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getElementaryCharacteristics",
        &PhysicalProblem::getElementaryCharacteristics, R"(
Return the elementary charateristics

Returns:
    ElementaryCharacteristicsPtr: a pointer to the elementary charateristics
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getDOFNumbering", &PhysicalProblem::getDOFNumbering, R"(
Return the DOF numbering

Returns:
    BaseDOFNumberingPtr: a pointer to the DOF numbering
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getBehaviourProperty", &PhysicalProblem::getBehaviourProperty, R"(
Return the behaviour properties

Returns:
    BehaviourPropertyPtr: a pointer to the behaviour properties
        )",
            ( py::arg( "self" ) ) );
    c1.def( "computeListOfLoads",
        &PhysicalProblem::computeListOfLoads, R"(
Build the list of loads from the added loads

Returns:
    Bool: True if success
        )",
              ( py::arg( "self" ) )   );
    c1.def( "computeDOFNumbering", &PhysicalProblem::computeDOFNumbering, R"(
Build DOF numbering from the model and loads

Returns:
    Bool: True if success
        )",
            ( py::arg( "self" ) ) );
    c1.def( "getListOfLoads",
        &PhysicalProblem::getListOfLoads, R"(
Return list of loads.

Returns:
    ListOfLoadsPtr: a pointer to list of loads
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getListOfDirichletBCs",
        &PhysicalProblem::getListOfDirichletBCs, R"(
Return list of DirichletBCs

Returns:
    ListDiriBC: a list of DirichletBC
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getListOfMechanicalLoadsReal",
        &PhysicalProblem::getListOfMechanicalLoadsReal, R"(
Return list of real mechanical loads

Returns:
    ListMecaLoadReal: a list of real mechanical loads
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getListOfMechanicalLoadsFunction",
        &PhysicalProblem::getListOfMechanicalLoadsFunction, R"(
Return list of Function mechanical loads

Returns:
    ListMecaLoadFunction: a list of Function mechanical loads
        )",
              ( py::arg( "self" ) )   );
#ifdef ASTER_HAVE_MPI
    c1.def( "getListOfParallelMechanicalLoadsReal",
        &PhysicalProblem::getListOfParallelMechanicalLoadsReal, R"(
Return list of real parallel mechanical loads

Returns:
    ListParaMecaLoadReal: a list of real parallel mechanical loads
        )",
              ( py::arg( "self" ) )   );
    c1.def( "getListOfParallelMechanicalLoadsFunction",
        &PhysicalProblem::getListOfParallelMechanicalLoadsFunction, R"(
Return list of function parallel mechanical loads

Returns:
    ListParaMecaLoadFunction: a list of function parallel mechanical loads
        )",
              ( py::arg( "self" ) )   );
#endif /* ASTER_HAVE_MPI */
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
#endif /* ASTER_HAVE_MPI */
    addThermalLoadToInterface( c1 );
};
