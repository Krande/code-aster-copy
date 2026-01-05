/**
 * @file PetscRedistribute.cxx
 * @brief Given a factored matrix, here is a petsc wrapping function
 *        to effciently compute from a sparse RHS matrixonly a subset
 *        of entries of solution matrix.
 * @author Nicolas Tardieu
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "Utilities/PetscApplyFactorOnSubBlocks.h"

#ifdef ASTER_HAVE_PETSC
/* Helper: build J = { rows having any nonzero in B } */
PetscErrorCode BuildJFromRHS( Mat B, IS *Jout ) {
    PetscInt N, nrhs;
    PetscInt *Jbuf = NULL, cap = 64, nJ = 0;

    PetscFunctionBegin;
    PetscCall( MatGetSize( B, &N, &nrhs ) );
    PetscCall( PetscMalloc1( cap, &Jbuf ) );

    for ( PetscInt i = 0; i < N; i++ ) {
        const PetscInt *cols = NULL;
        const PetscScalar *vals = NULL;
        PetscInt ncols = 0;
        PetscReal tol = 0;
        PetscCall( MatGetRow( B, i, &ncols, &cols, &vals ) );
        if ( ncols > 0 ) {

            // PetscReal nrmInf = 0.0;
            // for (PetscInt jj = 0; jj < ncols; ++jj) {
            // PetscReal a = PetscAbsScalar(vals[jj]);
            // if (a > nrmInf) nrmInf = a;
            // }

            // if (nrmInf > tol) {                /* or: if (nrm2 > tol) */
            if ( nJ == cap ) {
                cap = (PetscInt)( ( 3.0 / 2.0 ) * cap + 16 );
                PetscCall( PetscRealloc( sizeof( *Jbuf ) * cap, &Jbuf ) );
            }
            Jbuf[nJ++] = i; /* keep this row */
            // }
        }
        PetscCall( MatRestoreRow( B, i, &ncols, &cols, &vals ) );
    }
    PetscCall( ISCreateGeneral( PETSC_COMM_SELF, nJ, Jbuf, PETSC_OWN_POINTER, Jout ) );
    PetscCall( ISSortRemoveDups( *Jout ) );
    PetscFunctionReturn( PETSC_SUCCESS );
}

/* Build R using MatCreateSeqAIJWithArrays so only rows in I contain columns J */
PetscErrorCode BuildSelectedInversePattern( PetscInt N, IS Iall, const PetscInt *Iidx, PetscInt nI,
                                            IS Jall, const PetscInt *Jidx, PetscInt nJ,
                                            Mat *Rout ) {
    PetscInt *ia = NULL, *ja = NULL;
    PetscScalar *aa = NULL;
    PetscBool *inI = NULL;
    PetscInt pos = 0;

    PetscFunctionBegin;

    /* CSR arrays: ia length N+1, ja/aa length nI*nJ */
    PetscCall( PetscCalloc1( N, &inI ) );
    for ( PetscInt k = 0; k < nI; k++ )
        inI[Iidx[k]] = PETSC_TRUE;

    PetscCall( PetscMalloc1( N + 1, &ia ) );
    PetscCall( PetscMalloc1( (size_t)nI * (size_t)nJ, &ja ) );
    PetscCall( PetscMalloc1( (size_t)nI * (size_t)nJ, &aa ) );

    ia[0] = 0;
    for ( PetscInt r = 0; r < N; r++ ) {
        if ( inI[r] ) {
            /* copy J into this row */
            for ( PetscInt j = 0; j < nJ; j++ ) {
                ja[pos] = Jidx[j];
                aa[pos] = 1.0;
                pos++;
            }
        }
        ia[r + 1] = pos;
    }
    PetscCall( PetscFree( inI ) );
    /* pos should equal nI*nJ */
    PetscCheck( pos == nI * nJ, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                "nnz mismatch: pos=%" PetscInt_FMT " != nI*nJ=%" PetscInt_FMT, pos, nI * nJ );

    /* Create R with user-provided CSR arrays (sequential) */
    PetscCall( MatCreateSeqAIJWithArrays( PETSC_COMM_SELF, N, N, ia, ja, aa, Rout ) );

    PetscFunctionReturn( PETSC_SUCCESS );
}

/**
 * @brief Solve a multi-RHS system using a MUMPS factor and a sparse PETSc RHS.
 *
 * Given a MUMPS factor matrix and a sparse PETSc Mat as the right-hand side,
 * this routine solves for all columns. If an index set I is provided, it first
 * determines the set J of row indices where the RHS has nonzero entries, and
 * returns a sparse MATAIJ solution containing only the I×J block. Supplying I
 * can speed up the solve, but is most effective when the RHS column count is
 * comparable to or larger than |J|.
 */
PetscErrorCode MatSparseSolve_petsc( Mat FctMat, Mat RhsMat, IS ISet, Mat *SolMat, IS *JSet ) {
    PetscInt N, N2, nrowsB, nrhs;
    PetscLogDouble t0, dt, dtmax;
    IS Iall = NULL, Jall = NULL, AllCols = NULL;
    const PetscInt *Iidx = NULL, *Jidx = NULL;
    PetscInt nI = 0, nJ = 0;
    Mat R = NULL, Rt = NULL, Y = NULL, S = NULL;
    MPI_Comm ca, cb;
    PetscMPIInt sizeA, sizeB, rank;

    PetscFunctionBegin;

    PetscCall( PetscObjectGetComm( (PetscObject)FctMat, &ca ) );
    PetscCall( PetscObjectGetComm( (PetscObject)RhsMat, &cb ) );
    PetscCallMPI( MPI_Comm_size( ca, &sizeA ) );
    PetscCallMPI( MPI_Comm_size( cb, &sizeB ) );
    PetscCheck(
        sizeA == 1 && sizeB == 1, PETSC_COMM_SELF, PETSC_ERR_SUP,
        "Sequential-only implementation: FctMat and RhsMat must live on a size-1 communicator" );

    PetscCall( MatGetSize( FctMat, &N, &N2 ) );
    PetscCall( MatGetSize( RhsMat, &nrowsB, &nrhs ) );
    PetscCheck( N == N2, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "FctMat must be square" );
    PetscCheck( nrowsB == N, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
                "RhsMat rows %" PetscInt_FMT " != %" PetscInt_FMT, nrowsB, N );

    /* I: keep as given or use all rows */
    if ( ISet ) {
        PetscCall( ISDuplicate( ISet, &Iall ) );
        PetscCall( ISSortRemoveDups( Iall ) );
    } else {
        PetscCall( ISCreateStride( PETSC_COMM_SELF, N, 0, 1, &Iall ) );
    }
    PetscCall( ISGetLocalSize( Iall, &nI ) );

    /* -----------------------------------------------------------------------------------*/
    /* ------------------------------- GetInverse APPROACH -------------------------------*/
    /* -----------------------------------------------------------------------------------*/

    PetscCall( ISGetIndices( Iall, &Iidx ) );
    /* J: union of nonzero rows of B */
    PetscCall( BuildJFromRHS( RhsMat, &Jall ) );
    PetscCall( ISGetLocalSize( Jall, &nJ ) );
    PetscCall( ISGetIndices( Jall, &Jidx ) );

    /* Quick exits */
    if ( !nI || !nJ ) {
        PetscCall( MatCreateSeqAIJ( PETSC_COMM_SELF, nI, nrhs, 0, NULL, &S ) );
        PetscCall( MatAssemblyBegin( S, MAT_FINAL_ASSEMBLY ) );
        PetscCall( MatAssemblyEnd( S, MAT_FINAL_ASSEMBLY ) );
        *SolMat = S;
        PetscCall( ISRestoreIndices( Iall, &Iidx ) );
        PetscCall( ISRestoreIndices( Jall, &Jidx ) );
        PetscCall( ISDestroy( &Iall ) );
        PetscCall( ISDestroy( &Jall ) );
        PetscFunctionReturn( PETSC_SUCCESS );
    }

    /* R pattern via MatCreateSeqAIJWithArrays (CSR) */
    PetscCallMPI( MPI_Barrier( PETSC_COMM_WORLD ) );
    t0 = MPI_Wtime();
    PetscCall( BuildSelectedInversePattern( N, Iall, Iidx, nI, Jall, Jidx, nJ, &R ) );

    /* Fill R with selected inverse entries */
    PetscCall( MatTranspose( R, MAT_INPLACE_MATRIX, &R ) );
    PetscCall( MatCreateTranspose( R, &Rt ) );

    dt = MPI_Wtime() - t0;
    PetscCallMPI( MPI_Reduce( &dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD ) );
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    if ( !rank ) {
        PetscCall( PetscPrintf( PETSC_COMM_SELF,
                                "----------------------------------------------------------\n"
                                "------------------ Substructured Matrix ------------------\n"
                                "----------------- Construction Time Table ----------------\n"
                                "----------------------------------------------------------\n"
                                "  Pre Solve time: %.6f seconds\n",
                                (double)dtmax ) );
    }

    PetscCallMPI( MPI_Barrier( PETSC_COMM_WORLD ) );
    t0 = MPI_Wtime();
    PetscCall( MatMumpsGetInverse( FctMat, Rt ) );
    dt = MPI_Wtime() - t0;
    /* Timing*/
    dt = MPI_Wtime() - t0;
    PetscCallMPI( MPI_Reduce( &dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD ) );
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    if ( !rank ) {
        PetscCall(
            PetscPrintf( PETSC_COMM_SELF, "  Mat Solve time: %.6f seconds\n", (double)dtmax ) );
    }
    PetscCallMPI( MPI_Barrier( PETSC_COMM_WORLD ) );
    t0 = MPI_Wtime();
    /* Reset ICNTL(30) so later normal solves won’t require NRHS=N */
    PetscCall( MatMumpsSetIcntl( FctMat, 30, 0 ) );
    PetscCall( MatDestroy( &Rt ) );

    /* Y = R * B */
    PetscCall( MatTransposeMatMult( R, RhsMat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y ) );
    /* S = Y[I,:] */
    PetscCall( ISCreateStride( PETSC_COMM_SELF, nrhs, 0, 1, &AllCols ) );
    PetscCall( MatCreateSubMatrix( Y, Iall, AllCols, MAT_INITIAL_MATRIX, &S ) );
    PetscCall( ISRestoreIndices( Iall, &Iidx ) );
    PetscCall( ISRestoreIndices( Jall, &Jidx ) );
    /* Clean up */
    PetscCall( MatDestroy( &Y ) );
    PetscCall( MatDestroy( &R ) );
    PetscCall( ISDestroy( &Iall ) );
    PetscCall( ISDestroy( &Jall ) );
    PetscCall( ISDestroy( &AllCols ) );

    dt = MPI_Wtime() - t0;
    PetscCallMPI( MPI_Reduce( &dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD ) );
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );

    if ( !rank ) {
        PetscCall( PetscPrintf( PETSC_COMM_SELF, "  Post Solve time: %.6f s \n", (double)dtmax ) );
    }

    /* -----------------------------------------------------------------------------*/

    /* Building the final parallel matrix */
    if ( JSet == NULL ) {
        /* Return */
        *SolMat = S;
    } else {
        PetscInt Ng = 0;
        Mat A_diag = NULL, A = NULL;

        PetscCallMPI( MPI_Barrier( PETSC_COMM_WORLD ) );
        t0 = MPI_Wtime();
        PetscCallMPI( MPI_Allreduce( &nI, &Ng, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD ) );
        PetscCheck( JSet && *JSet, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "JSet must be provided" );
        PetscLogDouble t00 = MPI_Wtime();
        PetscCall( MatCreateSeqAIJ( PETSC_COMM_SELF, nI, nI, 1, NULL, &A_diag ) );
        for ( PetscInt i = 0; i < nI; ++i )
            PetscCall( MatSetValue( A_diag, i, i, 1.0, INSERT_VALUES ) );
        PetscCall( MatAssemblyBegin( A_diag, MAT_FINAL_ASSEMBLY ) );
        PetscCall( MatAssemblyEnd( A_diag, MAT_FINAL_ASSEMBLY ) );
        dt = MPI_Wtime() - t00;
        PetscCallMPI( MPI_Reduce( &dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD ) );
        if ( !rank )
            PetscCall( PetscPrintf( PETSC_COMM_SELF, "  A_diag time: %.6f s \n", (double)dtmax ) );
        const PetscInt *JsetIdx = NULL;
        PetscCall( ISGetIndices( *JSet, &JsetIdx ) );
        PetscCall( MatCreateMPIAIJWithSeqAIJ( PETSC_COMM_WORLD, Ng, Ng, A_diag, S,
                                              (PetscInt *)JsetIdx, &A ) );
        PetscCall( ISRestoreIndices( *JSet, &JsetIdx ) );
        PetscCall( MatSetOption( A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE ) );
        PetscCall( MatSetOption( A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE ) );
        PetscCall( MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ) );
        PetscCall( MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ) );
        *SolMat = A;
        dt = MPI_Wtime() - t0;
        PetscCallMPI( MPI_Barrier( PETSC_COMM_WORLD ) );
        PetscCallMPI( MPI_Reduce( &dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD ) );
        if ( !rank )
            PetscCall( PetscPrintf( PETSC_COMM_SELF,
                                    "  Assembling time: %.6f s \n"
                                    "----------------------------------------------------------\n",
                                    (double)dtmax ) );
    }
    PetscFunctionReturn( PETSC_SUCCESS );
}
#else
void MatSparseSolve_petsc() { std::cout << "PETSc library non available" << std::endl; };
#endif

/**
 * @brief This the python binding function to MatSparseSolve_petsc
 */
py::object applyFactorOnSubBlocks( py::object pyFctMat, py::object pyRhsMat, py::object pyISet,
                                   py::object pyJSet ) {
#ifdef ASTER_HAVE_PETSC4PY
    PetscErrorCode ierr;
    Mat new_mat;

    /* ------------------ Recasting Python objects to C objects ------------------*/

    /* Mats */
    Mat FctMat, RhsMat;
    struct PyPetscMatObject *pyx_fctMat = (struct PyPetscMatObject *)( pyFctMat.ptr() );
    FctMat = pyx_fctMat->mat;
    struct PyPetscMatObject *pyx_rhsMat = (struct PyPetscMatObject *)( pyRhsMat.ptr() );
    RhsMat = pyx_rhsMat->mat;

    /* IS */
    IS ISet, JSet;
    struct PyPetscISObject *pyx_iSet = (struct PyPetscISObject *)( pyISet.ptr() );
    ISet = pyx_iSet->iset;
    if ( pyJSet.is_none() ) {
        /* C function call */
        ierr = MatSparseSolve_petsc( FctMat, RhsMat, ISet, &new_mat );
    } else {
        struct PyPetscISObject *pyx_jSet = (struct PyPetscISObject *)( pyJSet.ptr() );
        JSet = pyx_jSet->iset;
        /* C function call */
        ierr = MatSparseSolve_petsc( FctMat, RhsMat, ISet, &new_mat, &JSet );
    }

    assert( ierr == 0 );
    if ( new_mat ) {
        // petsc4py shell for the new petsc4py matrix
        py::object petsc_matr = py::module_::import( "petsc4py.PETSc" ).attr( "Mat" )();
        // Recasting output Mat as a petsc4py.Mat
        struct PyPetscMatObject *pyx_mat = (struct PyPetscMatObject *)( petsc_matr.ptr() );
        pyx_mat->mat = new_mat;
        py::object result = py::reinterpret_steal< py::object >( (PyObject *)pyx_mat );
        petsc_matr.release();
        // Return the output petsc4py matrix
        return result;
    } else {
        return py::none();
    }
#else
    PyErr_SetString( PyExc_NotImplementedError, "petsc4py is not available" );
    throw py::error_already_set();
#endif
}
