/**
 * @file PhysicalProblemInterface.cxx
 * @brief Interface python de PhysicalProblem
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

#include "PythonBindings/PhysicalProblemInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/BaseDOFNumberingInterface.h"
#include "PythonBindings/LoadUtilities.h"

void exportPhysicalProblemToPython( py::module_ &mod ) {

    py::class_< PhysicalProblem, PhysicalProblemPtr > c1( mod, "PhysicalProblem" );
    // fake initFactoryPtr: not a DataStructure
    c1.def( py::init( &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr > ) );
    c1.def( py::init( &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr,
                                       ElementaryCharacteristicsPtr > ) );
    c1.def( "getModel", &PhysicalProblem::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )" );
    c1.def( "getMesh", &PhysicalProblem::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )" );
    c1.def( "getMaterialField", &PhysicalProblem::getMaterialField, R"(
Return the material field

Returns:
    MaterialFieldPtr: a pointer to the material field
        )" );
    c1.def( "getCodedMaterial", &PhysicalProblem::getCodedMaterial, R"(
Return the coded material

Returns:
    CodedMaterialPtr: a pointer to the coded material
        )" );
    c1.def( "getElementaryCharacteristics", &PhysicalProblem::getElementaryCharacteristics, R"(
Return the elementary charateristics

Returns:
    ElementaryCharacteristicsPtr: a pointer to the elementary charateristics
        )" );
    c1.def( "getDOFNumbering", &PhysicalProblem::getDOFNumbering, R"(
Return the DOF numbering

Returns:
    BaseDOFNumberingPtr: a pointer to the DOF numbering
        )" );
    c1.def( "setDOFNumbering", &PhysicalProblem::setDOFNumbering, R"(
Set the DOF numbering

Arguments:
    BaseDOFNumberingPtr: a pointer to the DOF numbering
        )",
            py::arg( "dofNume" ) );
    c1.def( "getBehaviourProperty", &PhysicalProblem::getBehaviourProperty, R"(
Return the behaviour properties

Returns:
    BehaviourPropertyPtr: a pointer to the behaviour properties
        )" );
    c1.def( "computeListOfLoads", &PhysicalProblem::computeListOfLoads, R"(
Build the list of loads from the added loads

Returns:
    Bool: True if success
        )" );
    c1.def( "computeDOFNumbering", &PhysicalProblem::computeDOFNumbering, R"(
Build DOF numbering from the model and loads

Returns:
    Bool: True if success
        )" );
    c1.def( "computeBehaviourProperty", &PhysicalProblem::computeBehaviourProperty,
            R"(
    Create constant fields on cells for behaviour (COMPOR, CARCRI and MULCOM)

    Arguments:
        COMPORTEMENT (list[dict]): keywords as provided to STAT_NON_LINE/COMPORTEMENT
        SIGM_INIT (str): "OUI" if there is an initial stress field
        INFO (int): level of verbosity, 1 to have description of behaviour or 0 to be quiet
            )",
            py::arg( "COMPORTEMENT" ), py::arg( "SIGM_INIT" ) = "NON", py::arg( "INFO" ) = 1 );
    c1.def( "getListOfLoads", &PhysicalProblem::getListOfLoads, R"(
Return list of loads.

Returns:
    ListOfLoadsPtr: a pointer to list of loads
        )" );
    c1.def( "getExternalStateVariables", &PhysicalProblem::getExternalStateVariables, R"(
    Get the field for external state variables

    Arguments:
        time [float] : time value to evaluate values

    Returns:
        FieldOnCellsRealPtr : external values
          )",
            py::arg( "time" ) );
    c1.def( "getReferenceExternalStateVariables",
            &PhysicalProblem::getReferenceExternalStateVariables, R"(
    Get the field of reference values for external state variables

    Returns:
        FieldOnCellsRealPtr : field of reference values
          )" );
    c1.def( "computeReferenceExternalStateVariables",
            &PhysicalProblem::computeReferenceExternalStateVariables, R"(
    Compute field for external state variables reference value

    Returns:
        FieldOnCells: field for external state variables reference values
            )" );
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
    addParallelThermalLoadToInterface( c1 );
#endif /* ASTER_HAVE_MPI */
    addThermalLoadToInterface( c1 );
    addAcousticLoadToInterface( c1 );
};
