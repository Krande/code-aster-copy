/**
 * @file PhysicalProblemInterface.cxx
 * @brief Interface python de PhysicalProblem
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

#include "PythonBindings/PhysicalProblemInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/BaseDOFNumberingInterface.h"
#include "PythonBindings/LoadUtilities.h"

void exportPhysicalProblemToPython( py::module_ &mod ) {

    py::class_< PhysicalProblem, PhysicalProblemPtr > c1( mod, "PhysicalProblem" );
    // fake initFactoryPtr: not a DataStructure
    c1.def( py::init( &initFactoryPtr< PhysicalProblem, BaseDOFNumberingPtr > ) );
    c1.def( py::init( &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr > ) );
    c1.def( py::init( &initFactoryPtr< PhysicalProblem, ModelPtr, MaterialFieldPtr,
                                       ElementaryCharacteristicsPtr > ) );
    c1.def( define_pickling< PhysicalProblem >() );

    c1.def( "getModel", &PhysicalProblem::getModel, R"(
Return the model

Returns:
    Model: a pointer to the model
        )" );
    c1.def( "getMesh", &PhysicalProblem::getMesh, R"(
Return the mesh

Returns:
    Mesh: a pointer to the mesh
        )" );
    c1.def( "getMaterialField", &PhysicalProblem::getMaterialField, R"(
Return the material field

Returns:
    MaterialField: a pointer to the material field
        )" );
    c1.def( "getCodedMaterial", &PhysicalProblem::getCodedMaterial, R"(
Return the coded material

Returns:
    CodedMaterial: a pointer to the coded material
        )" );
    c1.def( "getElementaryCharacteristics", &PhysicalProblem::getElementaryCharacteristics, R"(
Return the elementary charateristics

Returns:
    ElementaryCharacteristics: a pointer to the elementary charateristics
        )" );
    c1.def( "getDOFNumbering", &PhysicalProblem::getDOFNumbering, R"(
Return the DOF numbering

Returns:
    BaseDOFNumbering: a pointer to the DOF numbering
        )" );
    c1.def( "setDOFNumbering", &PhysicalProblem::setDOFNumbering, R"(
Set the DOF numbering

Arguments:
    dofNum (BaseDOFNumbering): a pointer to the DOF numbering
        )",
            py::arg( "dofNum" ) );
    c1.def( "getBehaviourProperty", &PhysicalProblem::getBehaviourProperty, R"(
Return the behaviour properties

Returns:
    BehaviourProperty: a pointer to the behaviour properties
        )" );
    c1.def( "computeListOfLoads", &PhysicalProblem::computeListOfLoads, R"(
Build the list of loads from the added loads

Arguments:
    command_name (str): It is possible to add a command name to add more checking (default: "")

Returns:
    Bool: True if success
        )",
            py::arg( "command_name" ) = std::string() );
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
    ListOfLoads: a pointer to list of loads
        )" );
    c1.def( "setListOfLoads", &PhysicalProblem::setListOfLoads, R"(
Set list of loads

Arguments:
    loads (ListOfLoads): a pointer to the list of loads
        )",
            py::arg( "loads" ) );
    c1.def( "getExternalStateVariables", &PhysicalProblem::getExternalStateVariables, R"(
    Get the field for external state variables

    Arguments:
        time [float] : time value to evaluate values

    Returns:
        FieldOnCellsReal : external values
          )",
            py::arg( "time" ) );
    c1.def( "getReferenceExternalStateVariables",
            &PhysicalProblem::getReferenceExternalStateVariables, R"(
    Get the field of reference values for external state variables

    Returns:
        FieldOnCellsReal : field of reference values
          )" );
    c1.def( "computeReferenceExternalStateVariables",
            &PhysicalProblem::computeReferenceExternalStateVariables, R"(
    Compute field for external state variables reference value

    Returns:
        FieldOnCells: field for external state variables reference values
            )" );
    // -----------------------------------------------------------------------------------------
    c1.def( "getDirichletBCDOFs", &PhysicalProblem::getDirichletBCDOFs, R"(
    Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

    Returns:
        tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
            size = neq + 1,
            tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
            tuple(neq) = number of DOFs eliminated.
        )" );
    // -----------------------------------------------------------------------------------------
    c1.def( "zeroDirichletBCDOFs", &PhysicalProblem::zeroDirichletBCDOFs, R"(
    Set in-place to zero the DOFs with DirichletBC (aka not assigned by Lagrange multipliers)

    Returns:
        field(FieldOnNodes): the modified field
        )" );
    c1.def( "isMechanical", &PhysicalProblem::isMechanical, R"(
            To know if the problem is mechanical or not

            Returns:
                bool: True - if the model is mechanical
            )" );
    c1.def( "isThermal", &PhysicalProblem::isThermal, R"(
            To know if the problem is thermal or not

            Returns:
                bool: True - if the model is thermal
            )" );
    c1.def( "isAcoustic", &PhysicalProblem::isAcoustic, R"(
            To know if the probleme is acoustic or not

            Returns:
                bool: True - if the model is acoustic
            )" );
    // -----------------------------------------------------------------------------------------
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
    c1.def( "setVirtualSlavCell", &PhysicalProblem::setVirtualSlavCell, R"(
        Set virtual cells from contact definition

        Arguments:
            virtualCell (FiniteElementDescriptor)): a pointer to the FED
                )",
            py::arg( "contact" ) );
    c1.def( "setVirtualCell", &PhysicalProblem::setVirtualCell, R"(
        Set virtual cells from contact pairing

        Arguments:
            virtualCell (FiniteElementDescriptor)): a pointer to the FED
                )",
            py::arg( "virtualCell" ) );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
    addParallelThermalLoadToInterface( c1 );
#endif /* ASTER_HAVE_MPI */
    addThermalLoadToInterface( c1 );
    addAcousticLoadToInterface( c1 );
};
