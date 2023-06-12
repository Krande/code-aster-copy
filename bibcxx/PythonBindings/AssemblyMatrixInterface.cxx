/**
 * @file AssemblyMatrixInterface.cxx
 * @brief Interface python de AssemblyMatrix
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/AssemblyMatrixInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/FieldOnNodesInterface.h"

void exportAssemblyMatrixToPython( py::module_ &mod ) {

    py::class_< AssemblyMatrixDisplacementReal, AssemblyMatrixDisplacementRealPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixDisplacementReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal, PhysicalProblemPtr > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal,
                                         const AssemblyMatrixDisplacementReal & > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixDisplacementReal::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixDisplacementReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixDisplacementReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "scale", &AssemblyMatrixDisplacementReal::scale, R"(
Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

Arguments:
    lvect (list[float]): List of the values.
    rvect (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixDisplacementReal::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "size", &AssemblyMatrixDisplacementReal::size, R"(
Get the size of the matrix

Arguments:
    local (bool) local or global size
        )",
              py::arg( "local" ) = true )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixDisplacementReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixDisplacementReal::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixDisplacementReal::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def( "__mul__", +[]( const AssemblyMatrixDisplacementReal &M, const FieldOnNodesReal &v ) {
            return M * v;
        } );

    py::class_< AssemblyMatrixEliminatedReal, AssemblyMatrixEliminatedRealPtr,
                AssemblyMatrixDisplacementReal >( mod, "AssemblyMatrixEliminatedReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixEliminatedReal > ) )
        .def( py::init( &initFactoryPtr< AssemblyMatrixEliminatedReal, std::string > ) );

    py::class_< AssemblyMatrixDisplacementComplex, AssemblyMatrixDisplacementComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixDisplacementComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixDisplacementComplex::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixDisplacementComplex]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixDisplacementComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixDisplacementComplex::transposeConjugate )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixDisplacementComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixDisplacementComplex::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementComplex::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixDisplacementComplex::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixDisplacementComplex::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixDisplacementComplex::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixDisplacementComplex::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( -py::self )
        .def( "__mul__", +[]( const AssemblyMatrixDisplacementComplex &M,
                              const FieldOnNodesComplex &v ) { return M * v; } );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureReal, AssemblyMatrixTemperatureRealPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixTemperatureReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal, PhysicalProblemPtr > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixTemperatureReal::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixTemperatureReal]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixTemperatureReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixTemperatureReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixTemperatureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "scale", &AssemblyMatrixTemperatureReal::scale, R"(
Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

Arguments:
    lvect (list[float]): List of the values.
    rvect (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixTemperatureReal::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "size", &AssemblyMatrixTemperatureReal::size, R"(
Get the size of the matrix

Arguments:
    local (bool) local or global size
        )",
              py::arg( "local" ) = true )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixTemperatureReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixTemperatureReal::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixTemperatureReal::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def( "__mul__", +[]( const AssemblyMatrixTemperatureReal &M, const FieldOnNodesReal &v ) {
            return M * v;
        } );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureComplex, AssemblyMatrixTemperatureComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixTemperatureComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixTemperatureComplex::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixTemperatureComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixTemperatureComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixTemperatureComplex::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixTemperatureComplex::transposeConjugate );

    py::class_< AssemblyMatrixPressureReal, AssemblyMatrixPressureRealPtr, BaseAssemblyMatrix >(
        mod, "AssemblyMatrixPressureReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixPressureReal::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixPressureReal::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureReal::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixPressureReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureReal::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixPressureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        .def( "copy", &AssemblyMatrixPressureReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def( "__mul__", +[]( const AssemblyMatrixPressureReal &M, const FieldOnNodesReal &v ) {
            return M * v;
        } );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixPressureComplex, AssemblyMatrixPressureComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixPressureComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "addElementaryMatrix", &AssemblyMatrixPressureComplex::addElementaryMatrix, R"(
Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem

Arguments:
    matr_elem [ElementaryMatrixPressureComplex]: elementary matrix to add
    coeff [float]: assembling factor (default = 1.0)
        )",
              py::arg( "matr_elem" ), py::arg( "coeff" ) = 1.0 )
        // -----------------------------------------------------------------------------------------
        .def( "clearElementaryMatrix", &AssemblyMatrixPressureComplex::clearElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "getMaterialField", &AssemblyMatrixPressureComplex::getMaterialField )
        // -----------------------------------------------------------------------------------------
        .def( "getNumberOfElementaryMatrix",
              &AssemblyMatrixPressureComplex::getNumberOfElementaryMatrix )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixPressureComplex::transposeConjugate )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixPressureComplex::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixPressureComplex::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixPressureComplex::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixPressureComplex::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def( "__mul__", +[]( const AssemblyMatrixPressureComplex &M,
                              const FieldOnNodesComplex &v ) { return M * v; } );
};
