/**
 * @file AssemblyMatrixInterface.cxx
 * @brief Interface python de AssemblyMatrix
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/AssemblyMatrixInterface.h"
#include <PythonBindings/factory.h>

void exportAssemblyMatrixToPython() {

    void ( AssemblyMatrixDisplacementReal::*c1 )( const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixDisplacementReal::addLoad;
    void ( AssemblyMatrixDisplacementReal::*c2 )( const DirichletBCPtr &currentLoad,
                                                            const FunctionPtr &func ) =
        &AssemblyMatrixDisplacementReal::addLoad;

    py::class_< AssemblyMatrixDisplacementReal, AssemblyMatrixDisplacementRealPtr,
                py::bases< DataStructure > >( "AssemblyMatrixDisplacementReal", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixDisplacementReal >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixDisplacementReal, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c1 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c2 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixDisplacementReal::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixDisplacementReal::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixDisplacementReal::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "getModel", &AssemblyMatrixDisplacementReal::getModel, R"(
Return the model.

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" )))
// -------------------------------------------------------------------------------------------------
        .def( "getMesh", &AssemblyMatrixDisplacementReal::getMesh, R"(
Return the mesh.

Returns:
    MeshPtr: a pointer to the mesh
        )",
              ( py::arg( "self" ) ) )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementReal::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "isEmpty", &AssemblyMatrixDisplacementReal::isEmpty, R"(
Test if the matrix is empty.

Returns:
    Bool: true if the matrix is empty
        )",
              ( py::arg( "self" )) )
// -------------------------------------------------------------------------------------------------
.def( "isFactorized", &AssemblyMatrixDisplacementReal::isFactorized, R"(
Test if the matrix is factorized.

Returns:
    Bool: true if the matrix is factorized
        )",
              ( py::arg( "self" )) )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementReal::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixDisplacementReal::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixDisplacementReal::setSolverName )
// -------------------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementReal::setValues, R"(
Erase the assembly matrix and set new values in it.
The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )")
// -------------------------------------------------------------------------------------------------
        .def( "hasDirichletEliminationDOFs",
            &AssemblyMatrixDisplacementReal::hasDirichletEliminationDOFs, R"(
Return True if matrix has some DOFs eliminated by Dirichlet boundaries conditions

Return:
    bool: True if matrix has some DOFs eliminated by Dirichlet boundaries conditions else False
        )")
// -------------------------------------------------------------------------------------------------
        .def( "getDirichletBCDOFs",
            &AssemblyMatrixDisplacementReal::getDirichletBCDOFs, R"(
Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

Return:
    tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions
            size = neq + 1

            tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0
            tuple(neq) = number of DOFs eliminated
        )")
// -------------------------------------------------------------------------------------------------
        .def( "getLagrangeScaling",
            &AssemblyMatrixDisplacementReal::getLagrangeScaling, R"(
Return the scaling used for Lagrange multipliers. It returns 1 if no Lagrange.

Return:
    double: scaling used for Lagrange multipliers. It returns 1 if no Lagrange
    are present.
        )")
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixDisplacementReal::*) () const>
                                                (&AssemblyMatrixDisplacementReal::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print",
            static_cast<void (AssemblyMatrixDisplacementReal::*) (const std::string) const>
                                                (&AssemblyMatrixDisplacementReal::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixDisplacementReal::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixDisplacementReal::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "transpose", &AssemblyMatrixDisplacementReal::transpose );


    void ( AssemblyMatrixDisplacementComplex::*c3 )(
        const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixDisplacementComplex::addLoad;
    void ( AssemblyMatrixDisplacementComplex::*c4 )( const DirichletBCPtr &currentLoad,
                                                             const FunctionPtr &func ) =
        &AssemblyMatrixDisplacementComplex::addLoad;

    py::class_< AssemblyMatrixDisplacementComplex, AssemblyMatrixDisplacementComplexPtr,
                py::bases< DataStructure > >( "AssemblyMatrixDisplacementComplex", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixDisplacementComplex >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixDisplacementComplex, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c3 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c4 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixDisplacementComplex::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixDisplacementComplex::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixDisplacementComplex::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "transpose", &AssemblyMatrixDisplacementComplex::transpose )
// -------------------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixDisplacementComplex::transposeConjugate )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementComplex::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementComplex::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixDisplacementComplex::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixDisplacementComplex::*) () const>
                                                (&AssemblyMatrixDisplacementComplex::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print",
            static_cast<void (AssemblyMatrixDisplacementComplex::*) (const std::string) const>
                                                (&AssemblyMatrixDisplacementComplex::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixDisplacementComplex::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixDisplacementComplex::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixDisplacementComplex::setSolverName );
// -------------------------------------------------------------------------------------------------


    void ( AssemblyMatrixTemperatureReal::*c5 )( const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixTemperatureReal::addLoad;
    void ( AssemblyMatrixTemperatureReal::*c6 )( const DirichletBCPtr &currentLoad,
                                                           const FunctionPtr &func ) =
        &AssemblyMatrixTemperatureReal::addLoad;

    py::class_< AssemblyMatrixTemperatureReal, AssemblyMatrixTemperatureRealPtr,
                py::bases< DataStructure > >( "AssemblyMatrixTemperatureReal", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixTemperatureReal >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixTemperatureReal, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c5 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c6 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixTemperatureReal::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixTemperatureReal::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixTemperatureReal::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "getModel", &AssemblyMatrixTemperatureReal::getModel, R"(
Return the model.

Returns:
    ModelPtr: a pointer to the model
        )",
              ( py::arg( "self" )))
// -------------------------------------------------------------------------------------------------
        .def( "getMesh", &AssemblyMatrixTemperatureReal::getMesh, R"(
Return the mesh.

Returns:
    MeshPtr: a pointer to the mesh
        )",
              ( py::arg( "self" ) ) )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureReal::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureReal::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixTemperatureReal::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixTemperatureReal::setSolverName )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixTemperatureReal::*) () const>
                                                (&AssemblyMatrixTemperatureReal::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print",
            static_cast<void (AssemblyMatrixTemperatureReal::*) (const std::string) const>
                                                (&AssemblyMatrixTemperatureReal::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixTemperatureReal::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixTemperatureReal::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixTemperatureReal::setValues, R"(
Erase the assembly matrix and set new values in it.
The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )");
// -------------------------------------------------------------------------------------------------


    void ( AssemblyMatrixTemperatureComplex::*c7 )( const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixTemperatureComplex::addLoad;
    void ( AssemblyMatrixTemperatureComplex::*c8 )( const DirichletBCPtr &currentLoad,
                                                            const FunctionPtr &func ) =
        &AssemblyMatrixTemperatureComplex::addLoad;

    py::class_< AssemblyMatrixTemperatureComplex, AssemblyMatrixTemperatureComplexPtr,
                py::bases< DataStructure > >( "AssemblyMatrixTemperatureComplex", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixTemperatureComplex >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixTemperatureComplex, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c7 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c8 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixTemperatureComplex::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixTemperatureComplex::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixTemperatureComplex::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureComplex::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureComplex::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixTemperatureComplex::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixTemperatureComplex::*) () const>
                                                (&AssemblyMatrixTemperatureComplex::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print",
            static_cast<void (AssemblyMatrixTemperatureComplex::*) (const std::string) const>
                                                (&AssemblyMatrixTemperatureComplex::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixTemperatureComplex::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixTemperatureComplex::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixTemperatureComplex::setSolverName );

    void ( AssemblyMatrixPressureReal::*c9 )( const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixPressureReal::addLoad;
    void ( AssemblyMatrixPressureReal::*c10 )( const DirichletBCPtr &currentLoad,
                                                         const FunctionPtr &func ) =
        &AssemblyMatrixPressureReal::addLoad;


    py::class_< AssemblyMatrixPressureReal, AssemblyMatrixPressureRealPtr,
                py::bases< DataStructure > >( "AssemblyMatrixPressureReal", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixPressureReal >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< AssemblyMatrixPressureReal, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c9 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c10 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixPressureReal::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixPressureReal::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixPressureReal::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureReal::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureReal::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixPressureReal::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixPressureReal::setSolverName )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixPressureReal::*) () const>
                                                (&AssemblyMatrixPressureReal::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixPressureReal::*) (const std::string) const>
                                                (&AssemblyMatrixPressureReal::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixPressureReal::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixPressureReal::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixPressureReal::setValues, R"(
Erase the assembly matrix and set new values in it.
The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )");
// -------------------------------------------------------------------------------------------------


    void ( AssemblyMatrixPressureComplex::*c11 )( const DirichletBCPtr &currentLoad ) =
        &AssemblyMatrixPressureComplex::addLoad;
    void ( AssemblyMatrixPressureComplex::*c12 )( const DirichletBCPtr &currentLoad,
                                                          const FunctionPtr &func ) =
        &AssemblyMatrixPressureComplex::addLoad;

    py::class_< AssemblyMatrixPressureComplex, AssemblyMatrixPressureComplexPtr,
                py::bases< DataStructure > >( "AssemblyMatrixPressureComplex", py::no_init )
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< AssemblyMatrixPressureComplex >))
// -------------------------------------------------------------------------------------------------
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< AssemblyMatrixPressureComplex, std::string >))
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c11 )
// -------------------------------------------------------------------------------------------------
        .def( "addDirichletBC", c12 )
// -------------------------------------------------------------------------------------------------
        .def( "appendElementaryMatrix",
              &AssemblyMatrixPressureComplex::appendElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "build", &AssemblyMatrixPressureComplex::build )
// -------------------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &AssemblyMatrixPressureComplex::getDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureComplex::getMaterialField )
// -------------------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureComplex::getNumberOfElementaryMatrix )
// -------------------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &AssemblyMatrixPressureComplex::setDOFNumbering )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixPressureComplex::*) () const>
                                                (&AssemblyMatrixPressureComplex::print), R"(
Print the matrix in code_aster format (with information on the DOF).
        )",
              ( py::arg( "self" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print",
            static_cast<void (AssemblyMatrixPressureComplex::*) (const std::string) const>
                                                (&AssemblyMatrixPressureComplex::print), R"(
Print the matrix in the given format format.

Arguments:
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::arg( "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "print", static_cast<void (AssemblyMatrixPressureComplex::*) (const ASTERINTEGER,
                        const std::string) const> (&AssemblyMatrixPressureComplex::print), R"(
Print the matrix in the given format format and in the given logical unit.

Arguments:
    unit (int) : logical unit to print
    format (string): 'ASTER' or 'MATLAB'

        )",
              ( py::arg( "self" ), py::args( "unit", "format" ))  )
// -------------------------------------------------------------------------------------------------
        .def( "setSolverName", &AssemblyMatrixPressureComplex::setSolverName );
// -------------------------------------------------------------------------------------------------
};
