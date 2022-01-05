/**
 * @file AssemblyMatrixInterface.cxx
 * @brief Interface python de AssemblyMatrix
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/BaseAssemblyMatrixInterface.h"
#include <PythonBindings/factory.h>

void exportBaseAssemblyMatrixToPython() {

    void ( BaseAssemblyMatrix::*c1 )( const DirichletBCPtr &currentLoad ) =
        &BaseAssemblyMatrix::addLoad;
    void ( BaseAssemblyMatrix::*c2 )( const DirichletBCPtr &currentLoad,
                                                  const FunctionPtr &func ) =
        &BaseAssemblyMatrix::addLoad;

    py::class_< BaseAssemblyMatrix, BaseAssemblyMatrixPtr,
                py::bases< DataStructure > >( "BaseAssemblyMatrix", py::no_init )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< BaseAssemblyMatrix, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< BaseAssemblyMatrix, std::string, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< BaseAssemblyMatrix, PhysicalProblemPtr, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addDirichletBC", c1 )
        // -----------------------------------------------------------------------------------------
        .def( "addDirichletBC", c2 )
        // -----------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &BaseAssemblyMatrix::getDOFNumbering )
        // -----------------------------------------------------------------------------------------
        .def( "getModel", &BaseAssemblyMatrix::getModel, R"(
Return the model.

Returns:
    Model: a pointer to the model
        )",
              ( py::arg( "self" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "getMesh", &BaseAssemblyMatrix::getMesh, R"(
Return the mesh.

Returns:
    Mesh: a pointer to the mesh
        )",
              ( py::arg( "self" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "getListOfLoads", &BaseAssemblyMatrix::getListOfLoads, R"(
Return the list of loads.

Returns:
    ListOfLoads: a pointer to the list of loads
        )",
              ( py::arg( "self" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "setListOfLoads", &BaseAssemblyMatrix::setListOfLoads, R"(
Set the list of loads.

Arguments:
    ListOfLoads: a pointer to the list of loads to set
        )",
              ( py::arg( "self" ), py::arg( "load" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "isEmpty", &BaseAssemblyMatrix::isEmpty, R"(
Test if the matrix is empty.

Returns:
    bool: *True* if the matrix is empty.
        )",
              ( py::arg( "self" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "isFactorized",
              static_cast< bool ( BaseAssemblyMatrix::* )() const >(
                  &BaseAssemblyMatrix::isFactorized ),
              R"(
Test if the matrix is factorized.

Returns:
    bool: *True* if the matrix is factorized.
        )",
              ( py::arg( "self" ) ) )
        .def( "isFactorized",
              static_cast< void ( BaseAssemblyMatrix::* )( const bool & ) >(
                  &BaseAssemblyMatrix::isFactorized ),
              R"(
Tell if the matrix is factorized.

Arguments:
    bool: *True* if the matrix is factorized.
        )",
              ( py::arg( "self" ), py::arg( "facto" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &BaseAssemblyMatrix::setDOFNumbering )
        // -----------------------------------------------------------------------------------------
        .def( "setSolverName", &BaseAssemblyMatrix::setSolverName )
        // -----------------------------------------------------------------------------------------
        .def( "hasDirichletEliminationDOFs",
              &BaseAssemblyMatrix::hasDirichletEliminationDOFs, R"(
Tell if matrix has some DOFs eliminated by Dirichlet boundaries conditions.

Returns:
    bool: *True* if matrix has some DOFs eliminated by Dirichlet boundaries conditions else *False*
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getDirichletBCDOFs", &BaseAssemblyMatrix::getDirichletBCDOFs, R"(
Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

Returns:
    tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
        size = neq + 1,
        tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
        tuple(neq) = number of DOFs eliminated.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getLagrangeScaling", &BaseAssemblyMatrix::getLagrangeScaling, R"(
Return the scaling used for Lagrange multipliers. It returns 1 if no Lagrange.

Returns:
    float: scaling used for Lagrange multipliers. It returns 1 if no Lagrange
    are present.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "print",
              static_cast< void ( BaseAssemblyMatrix::* )() const >(
                  &BaseAssemblyMatrix::print ),
              R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "print",
              static_cast< void ( BaseAssemblyMatrix::* )( const std::string ) const >(
                  &BaseAssemblyMatrix::print ),
              R"(
Print the matrix in the given format format.

Arguments:
    format (str): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "print",
              static_cast< void ( BaseAssemblyMatrix::* )( const ASTERINTEGER,
                                                                       const std::string ) const >(
                  &BaseAssemblyMatrix::print ),
              R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int): logical unit to print
    format (str): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ) ) )
        // -----------------------------------------------------------------------------------------
        .def( "transpose", &BaseAssemblyMatrix::transpose );
};
