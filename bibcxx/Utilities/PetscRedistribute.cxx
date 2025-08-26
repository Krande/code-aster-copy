/**
 * @file PetscRedistribute.cxx
 * @brief Inspired from petsc PCTelescope, here is some wrapping functions
 *        that helps repartioning a matrix on a subcommunicator.
 * @author Nicolas Tardieu
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
#include "Utilities/PetscRedistribute.h"

/**
 * @brief Based on a subcommunicator and the global number of rows,
 *        this function returns the rows index set of each proc.
 */
#ifdef ASTER_HAVE_PETSC4PY
PetscErrorCode GetRowsIS( MPI_Comm comm, PetscSubcomm psubcomm, PetscInt nr, IS *isin ) {

    // Create the Rows index sets on each proc
    if ( psubcomm->color == 0 ) {
        Vec xred;
        PetscInt st, ed;
        MPI_Comm subcomm;
        subcomm = PetscSubcommChild( psubcomm );
        // Create a redistributed model vector on the subcommunicator
        PetscCall( VecCreate( subcomm, &xred ) );
        PetscCall( VecSetType( xred, VECMPI ) );
        PetscCall( VecSetFromOptions( xred ) );
        PetscCall( VecSetSizes( xred, PETSC_DECIDE, nr ) );
        PetscCall( VecGetOwnershipRange( xred, &st, &ed ) );
        PetscCall( ISCreateStride( comm, ed - st, st, 1, isin ) );
    } else {
        PetscCall( ISCreateStride( comm, 0, 0, 1, isin ) );
    }

    PetscFunctionReturn( 0 );
};
#endif

/**
 * @brief Given a PETSc distributed matrix on an MPI communicator, this
 *        function gives back a redistributed matrix on a subcommunicator.
 */

py::object redistributePetscMat( py::object pMat, int subCommSize ) {
#ifdef ASTER_HAVE_PETSC4PY

    PetscErrorCode ierr;
    Mat new_mat;
    struct PyPetscMatObject *pyx_mat = (struct PyPetscMatObject *)( pMat.ptr() );
    Mat mat = pyx_mat->mat;
    ierr = redistribute_petsc( mat, subCommSize, &new_mat );
    assert( ierr == 0 );

    if ( new_mat ) {
        // petsc4py shell for the new petsc4py matrix
        py::object petsc_matr = py::module_::import( "petsc4py.PETSc" ).attr( "Mat" )();
        // Recasting redistributed Mat as a petsc4py.Mat
        pyx_mat = (struct PyPetscMatObject *)( petsc_matr.ptr() );
        pyx_mat->mat = new_mat;
        py::object result = py::reinterpret_steal< py::object >( (PyObject *)pyx_mat );
        petsc_matr.release();
        // Return the redistributed petsc4py matrix
        return result;
    } else {
        return py::none();
    }
#else
    PyErr_SetString( PyExc_NotImplementedError, "petsc4py is not available" );
    throw py::error_already_set();
#endif
}

#ifdef ASTER_HAVE_PETSC
static PetscErrorCode redistribute_petsc( Mat mat, int subCommSize, Mat *new_mat ) {

    // Variables
    MPI_Comm comm, subcomm;
    PetscSubcomm psubcomm;
    PetscInt nr, nc, bs;
    IS isrow, iscol;
    Mat Blocal, *_Blocal;
    // Getting COMM of the matrix to redistribute
    PetscObjectGetComm( (PetscObject)mat, &comm );
    // Defining the subcommunicator
    PetscCall( PetscSubcommCreate( comm, &psubcomm ) );
    PetscCall( PetscSubcommSetNumber( psubcomm, subCommSize ) );
    PetscCall( PetscSubcommSetType( psubcomm, PETSC_SUBCOMM_CONTIGUOUS ) );
    subcomm = PetscSubcommChild( psubcomm );
    // Global sizes
    PetscCall( MatGetSize( mat, &nr, &nc ) );
    // Rows index set
    PetscCall( GetRowsIS( comm, psubcomm, nr, &isrow ) );
    // Columns index set
    PetscCall( ISCreateStride( PETSC_COMM_SELF, nc, 0, 1, &iscol ) );
    // Each proc extracts its sequential part
    PetscCall( ISSetIdentity( iscol ) );
    PetscCall( MatGetBlockSizes( mat, NULL, &bs ) );
    PetscCall( ISSetBlockSize( iscol, bs ) );
    PetscCall( MatSetOption( mat, MAT_SUBMAT_SINGLEIS, PETSC_TRUE ) );
    PetscCall( MatCreateSubMatrices( mat, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, &_Blocal ) );
    Blocal = *_Blocal;
    PetscCall( PetscFree( _Blocal ) );
    // Creating the redistributed matrix
    if ( psubcomm->color == 0 ) {
        PetscInt mm;
        PetscCall( MatGetSize( Blocal, &mm, NULL ) );
        PetscCall(
            MatCreateMPIMatConcatenateSeqMat( subcomm, Blocal, mm, MAT_INITIAL_MATRIX, new_mat ) );
        // PetscCall(MatView(*new_mat, PETSC_VIEWER_STDOUT_(subcomm)));
    } else {
        ( *new_mat ) = NULL;
    }
    PetscCall( ISDestroy( &iscol ) );
    PetscCall( MatDestroy( &Blocal ) );

    PetscFunctionReturn( 0 );
}
#else
void redistribute_petsc() { std::cout << "PETSc library non available" << std::endl; };

#endif
