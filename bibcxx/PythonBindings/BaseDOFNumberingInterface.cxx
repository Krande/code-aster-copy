/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/BaseDOFNumberingInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/LoadUtilities.h"

void exportBaseDOFNumberingToPython( py::module_ &mod ) {

    py::class_< FieldOnNodesDescription, FieldOnNodesDescriptionPtr, DataStructure >(
        mod, "FieldOnNodesDescription" )
        .def( py::init( &initFactoryPtr< FieldOnNodesDescription > ) )
        .def( py::init( &initFactoryPtr< FieldOnNodesDescription, std::string > ) );

    py::class_< BaseDOFNumbering, BaseDOFNumbering::BaseDOFNumberingPtr, DataStructure > c1(
        mod, "BaseDOFNumbering" );
    // fake initFactoryPtr: created by subclasses
    // fake initFactoryPtr: created by subclasses
    c1.def( "addFiniteElementDescriptor", &BaseDOFNumbering::addFiniteElementDescriptor );
    c1.def( "computeNumbering", &BaseDOFNumbering::computeNumbering );
    c1.def( "getDescription", &BaseDOFNumbering::getDescription );
    c1.def( "getFiniteElementDescriptors", &BaseDOFNumbering::getFiniteElementDescriptors );
    c1.def( "getListOfLoads", &BaseDOFNumbering::getListOfLoads );
    c1.def( "setEmpty", &BaseDOFNumbering::setEmpty );
    c1.def( "getPhysicalQuantity", &BaseDOFNumbering::getPhysicalQuantity, R"(
Returns the name of the physical quantity that is numbered.

Returns:
    str: physical quantity name.
        )" );

    c1.def( "isParallel", &BaseDOFNumbering::isParallel, R"(
The numbering is distributed across MPI processes for High Performance Computing.

Returns:
    bool: *True* if used, *False* otherwise.
        )" );
    c1.def( "hasDirichletBC", &BaseDOFNumbering::hasDirichletBC, R"(
The list of loads used to build numbering contains Dirichlet BC.

Returns:
    bool: *True* if Dirichlet BC are present, *False* otherwise.
        )" );
    c1.def( "setElementaryMatrix", py::overload_cast< const ElementaryMatrixDisplacementRealPtr & >(
                                       &BaseDOFNumbering::setElementaryMatrix ) );
    c1.def( "setElementaryMatrix",
            py::overload_cast< const ElementaryMatrixDisplacementComplexPtr & >(
                &BaseDOFNumbering::setElementaryMatrix ) );
    c1.def( "setElementaryMatrix", py::overload_cast< const ElementaryMatrixTemperatureRealPtr & >(
                                       &BaseDOFNumbering::setElementaryMatrix ) );
    c1.def( "setElementaryMatrix", py::overload_cast< const ElementaryMatrixPressureComplexPtr & >(
                                       &BaseDOFNumbering::setElementaryMatrix ) );
    c1.def( "getModel", &BaseDOFNumbering::getModel, R"(
Return the model

Returns:
    ModelPtr: a pointer to the model
        )" );
    c1.def( "setModel", &BaseDOFNumbering::setModel );
    c1.def( "getMesh", &BaseDOFNumbering::getMesh, R"(
Return the mesh

Returns:
    MeshPtr: a pointer to the mesh
        )" );
    c1.def( "getDirichletBCDOFs", &BaseDOFNumbering::getDirichletBCDOFs, R"(
Return a vector which describes DOFs that are imposed by Dirichlet BC.

The vector has a size equals to the number of DOFs. For each dof, the value is equal to one
if Dirichel BC is imposed to this DOF else zero

Be carefull all Dirichlet BC have to be added before to call this function.

Returns:
    tuple(int): a list with dirichlet imposition.
        )" );
    addDirichletBCToInterface( c1 );
    addMechanicalLoadToInterface( c1 );
    addThermalLoadToInterface( c1 );
    addAcousticLoadToInterface( c1 );
#ifdef ASTER_HAVE_MPI
    addParallelMechanicalLoadToInterface( c1 );
#endif
};
