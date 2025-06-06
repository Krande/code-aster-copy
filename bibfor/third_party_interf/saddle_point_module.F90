! --------------------------------------------------------------------
! Copyright (C) 2016 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

module saddle_point_module
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
! person_in_charge: natacha.bereux at edf.fr
!
    use aster_petsc_module
    use saddle_point_context_type
!
    implicit none
!
    private
#include "asterc/asmpi_comm.h"
#include "asterfort/assert.h"

!
    public ::  convert_mat_to_saddle_point, update_double_lagrange, convert_rhs_to_saddle_point
!
#ifdef ASTER_HAVE_PETSC
    PetscErrorCode :: ierr
!
!
contains
!
!
! On entry a_mat is the PETsc matrix associated
! to Aster matrix matrasse. It contains a Double Lagrange system
! ( K C^T        C^T     )
! ( C -alpha I   alpha I )
! ( C  alpha I  -alpha I )
!
! On output, the matrix contains the following Simple Lagrange system
!
! a_mat = ( K C^T 0 )
!         ( C 0   0 )
!         ( 0 0   I )
! This matrix is defined as a MATSHELL matrix.
! K and C are stored in a context object.
!
    subroutine convert_mat_to_saddle_point(matrasse, a_mat)
        !
        use saddle_point_data_module, only: sp_context
        !
        ! Dummy arguments
        !
        character(len=19), intent(in) :: matrasse
        Mat, intent(inout)            :: a_mat
        !
        ! Local variables
        !
        ! Global size of the matrix
        PetscInt :: mg, ng
        ! Local size of the matrix
        PetscInt :: m, n
        mpi_int :: mpicomm
        PetscMPIInt :: comm
        !
        ! Récupération du communicateur MPI
        call asmpi_comm('GET', mpicomm)
        comm = mpicomm
        ! Init Saddle Point Context
        sp_context = new_saddle_point_context(matrasse, distributed_data, &
                                              a_mat)
        !
        ! The old matrix and the new one  shall have the same
        ! (global and local) sizes
        ! Get sizes
        call MatGetSize(a_mat, mg, ng, ierr)
        ASSERT(ierr == 0)
        ASSERT(mg == ng)
        call MatGetLocalSize(a_mat, m, n, ierr)
        ASSERT(ierr == 0)
        ASSERT(m == n)
        ! Destroy the old matrix
        call MatDestroy(a_mat, ierr)
        ASSERT(ierr == 0)
        !
        ! Create MatShell matrix
        ! TODO  utiliser l'objet  fortran context
        ! au lieu de PETSC_NULL_INTEGER
        call MatCreateShell(comm, m, m, mg, mg, PETSC_NULL_INTEGER, a_mat, ierr)
        ASSERT(ierr == 0)
        ! Define matrix-vector product operation
        call MatShellSetOperation(a_mat, MATOP_MULT, saddle_point_matmult, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine convert_mat_to_saddle_point
!
! On entry b is the right-hand side of a Double Lagrange system.
! On exit, it may be used with a Simple Lagrange System.
!
! The routine zeroes the section of b corresponding to Lagrange 2 multipliers
!
    subroutine convert_rhs_to_saddle_point(b)
        !
        use saddle_point_data_module, only: sp_context
        !
        ! Dummy arguments
        !
        Vec, intent(inout) :: b
        !
        ! Set to zero the Vector x3
        !
        call VecSet(sp_context%x3, 0.d0, ierr)
        ASSERT(ierr == 0)
        ! And load the values of x3 to vector b
        call VecScatterBegin(sp_context%scatter_to_lag2, sp_context%x3, b, &
          &     INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag2, sp_context%x3, b, &
         &         INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine convert_rhs_to_saddle_point
!
! MatMult operation used for the matshell matrix representing
! the simple Lagrange system matrix.
!
    subroutine saddle_point_matmult(a, x, y, perr)
        !
        use saddle_point_data_module, only: sp_context
        !
        ! Dummy arguments
        !
        Mat, intent(in)    :: a
        Vec, intent(in)    :: x
        Vec, intent(inout) :: y
        PetscErrorCode, intent(out) :: perr
        !
        ! Load values of  x1, x2 from x : use SCATTER_REVERSE
        call VecScatterBegin(sp_context%scatter_to_phys, x, sp_context%x1, &
        &     INSERT_VALUES, SCATTER_REVERSE, perr)
        ASSERT(perr == 0)
        call VecScatterEnd(sp_context%scatter_to_phys, x, sp_context%x1, &
            &    INSERT_VALUES, SCATTER_REVERSE, perr)
        ASSERT(perr == 0)
        !
        call VecScatterBegin(sp_context%scatter_to_lag1, x, sp_context%x2, &
            &    INSERT_VALUES, SCATTER_REVERSE, perr)
        ASSERT(perr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag1, x, sp_context%x2, &
            &    INSERT_VALUES, SCATTER_REVERSE, perr)
        ASSERT(perr == 0)
        !
        ! Do the Mat Mult :
        ! y1 = k_mat * x1 + c_mat^T * x2
        ! y2 = c_mat * x2
        !
        call MatMult(sp_context%k_mat, sp_context%x1, sp_context%xtmp, perr)
        ASSERT(perr == 0)
        call MatMultTransposeAdd(sp_context%c_mat, sp_context%x2, &
           &    sp_context%xtmp, sp_context%y1, perr)
        ASSERT(perr == 0)
        call MatMult(sp_context%c_mat, sp_context%x1, sp_context%y2, perr)
        ASSERT(perr == 0)
        !
        ! Load result to y : this time, use SCATTER_FORWARD
        !
        call VecScatterBegin(sp_context%scatter_to_phys, sp_context%y1, y, &
           &      INSERT_VALUES, SCATTER_FORWARD, perr)
        ASSERT(perr == 0)
        call VecScatterEnd(sp_context%scatter_to_phys, sp_context%y1, y, &
           &      INSERT_VALUES, SCATTER_FORWARD, perr)
        ASSERT(perr == 0)
        !
        call VecScatterBegin(sp_context%scatter_to_lag1, sp_context%y2, y, &
           &      INSERT_VALUES, SCATTER_FORWARD, perr)
        ASSERT(perr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag1, sp_context%y2, y, &
           &      INSERT_VALUES, SCATTER_FORWARD, perr)
        ASSERT(perr == 0)
        !
        perr = 0
        !
    end subroutine saddle_point_matmult
!
!
! On entry x contains:
! - the values of physical ddls
! - the values of (simple) Lagrange multipliers, stored
!   in x(islag1)
! On output:
!  - the values of physical ddls are unchanged
!  - the values of double Lagrange multipliers are computed
!    and returned in x(islag1) and x(islag2)
!
    subroutine update_double_lagrange(x)
        !
        use saddle_point_data_module, only: sp_context
        !
        ! Dummy arguments
        !
        Vec, intent(inout)    :: x
        !
        ! Local variables
        !
        PetscScalar :: alpha
        !
        ! Load values of x2 (containing simple lagrange value ) from x : use SCATTER_REVERSE
        !
        call VecScatterBegin(sp_context%scatter_to_lag1, x, sp_context%x2, &
             &     INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag1, x, sp_context%x2, &
             &    INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        !
        ! l1 = l2 = l/2,
        ! Scale lagrange multipliers ( x 0.5)
        alpha = 0.5d0
        call VecScale(sp_context%x2, alpha, ierr)
        ASSERT(ierr == 0)
        !
        ! Load result to lag1 and lag2 sections of x : this time, use SCATTER_FORWARD
        !
        call VecScatterBegin(sp_context%scatter_to_lag1, sp_context%x2, x, &
             &     INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag1, sp_context%x2, x, &
             &    INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterBegin(sp_context%scatter_to_lag2, sp_context%x3, x, &
             &    INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_context%scatter_to_lag2, sp_context%x3, x, &
             &    INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine update_double_lagrange
!
#else
contains
    subroutine convert_mat_to_saddle_point(matrasse, a_mat)
        character(len=19), intent(in) :: matrasse
        integer(kind=8), intent(in) :: a_mat
        character(len=1) :: kdummy
        integer(kind=8) :: idummy
        kdummy = matrasse(1:1)
        idummy = a_mat
    end subroutine convert_mat_to_saddle_point
!
    subroutine convert_rhs_to_saddle_point(b)
        integer(kind=8), intent(in) :: b
        integer(kind=8) :: idummy
        idummy = b
    end subroutine convert_rhs_to_saddle_point
!
    subroutine update_double_lagrange(x)
        integer(kind=8), intent(in)  :: x
        integer(kind=8) :: idummy
        idummy = x
    end subroutine update_double_lagrange
#endif
end module saddle_point_module
