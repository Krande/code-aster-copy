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

!
! A saddle_point_context object is a container used to manage
! a saddle point linear system
! ( k_mat c_mat^T ) (x_1) = (b_1)
! ( c_mat 0       ) (x_2)   (b_2)
! For simplicity the saddle point linear system is embedded
! in the larger system
! ( k_mat c_mat^T 0 ) (x_1) = (b_1)
! ( c_mat 0       0 ) (x_2)   (b_2)
! ( 0     0       Id) (x_3)   (b_3)
!
! It contains
! - Index Sets necessary to extract data from the global (double Lagrange) Aster system
! - matrix data ( k_mat, c_mat )
! - vector workspace ( x_1, x_2, b_1, b_2 )
!
module augmented_lagrangian_context_type
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use matrasse_module
    use saddle_point_context_type
!
    implicit none
!
    private
#include "asterf.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/assert.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/utmess.h"
!
!
    type, public :: augm_lagr_ctxt
        !
#ifdef ASTER_HAVE_PETSC
        type(saddlepoint_ctxt), pointer :: sp_ctxt => null()
        !
        ! Preconditioner data section
        ! ===========================
        ! Matrix used to build block (1,1) of the preconditioner
        Mat :: m_mat
        ! Scaling coefficient
        PetscReal :: gamma
        ! Preconditioner for the block of physical dofs
        PC :: pcphy
        !
#else
        integer(kind=8) :: idummy
#endif
    end type augm_lagr_ctxt
!
    public :: new_augmented_lagrangian_context, free_augm_lagrangian_context
!
#ifdef ASTER_HAVE_PETSC
    PetscErrorCode :: ierr
!
contains
!
!
    function new_augmented_lagrangian_context(sp_ctxt) result(ctxt)
        !
        type(saddlepoint_ctxt), target         :: sp_ctxt
        type(augm_lagr_ctxt)         :: ctxt
        !
        ! Local variables
        !
        ctxt%sp_ctxt => sp_ctxt
        !
        ! Init data section
        call set_precond_data(ctxt)
        !
    end function new_augmented_lagrangian_context
!
!
! Init data needed by an augmented lagrangian preconditioner
!
    subroutine set_precond_data(ctxt)
        !
        ! Dummy arguments
        !
        type(augm_lagr_ctxt), intent(inout) :: ctxt
        !
        !
        ! Local variables
        !
        mpi_int :: mpicomm
        mpi_int :: rang, nbproc
        PetscInt :: niremp
        PetscScalar :: fillp
        PetscReal :: aster_petsc_real
        type(saddlepoint_ctxt), pointer :: sp_ctxt => null()
        integer(kind=8), parameter :: icc_pre = 0, mumps_pre = 1
        integer(kind=8) :: pre_type
        Mat :: f
        !
        pre_type = mumps_pre
        ! TODO déterminer les bonnes valeurs
        niremp = 1
        fillp = 1.0d0
        !
        ! Notation
        sp_ctxt => ctxt%sp_ctxt
        !
        ! Récupération du communicateur MPI
        call asmpi_comm('GET', mpicomm)
        !
        ctxt%gamma = sp_ctxt%alpha
        !
        ! Compute block of physical dofs m_mat
        !
        ! m_mat = c^T c

        aster_petsc_real = PETSC_DEFAULT_REAL
        call MatTransposeMatMult(sp_ctxt%c_mat, sp_ctxt%c_mat,                  &
       &   MAT_INITIAL_MATRIX, aster_petsc_real, ctxt%m_mat, ierr)
        ASSERT(ierr == 0)
        ! m_mat <- k_mat + gamma*m_mat
        call MatAYPX(ctxt%m_mat, ctxt%gamma, sp_ctxt%k_mat, DIFFERENT_NONZERO_PATTERN, ierr)
        ASSERT(ierr == 0)
        call MatSetOption(ctxt%m_mat, MAT_SPD, PETSC_TRUE, ierr)
        ASSERT(ierr == 0)
        !
        ! Define preconditionner pcphy, based on the incomplete factorization of m_mat
        !
        call PCCreate(mpicomm, ctxt%pcphy, ierr)
        ASSERT(ierr == 0)
        call PCSetOperators(ctxt%pcphy, ctxt%m_mat, ctxt%m_mat, ierr)
        ASSERT(ierr == 0)
        !
        if (pre_type == icc_pre) then
            ! Attention, la factorisation ICC ne fonctionne qu'en séquentiel
            call asmpi_info(rank=rang, size=nbproc)
            if (nbproc > 1) then
                call utmess('F', 'PETSC_20')
            end if
            call PCSetType(ctxt%pcphy, PCICC, ierr)
            ASSERT(ierr == 0)
            call PCFactorSetLevels(ctxt%pcphy, niremp, ierr)
            ASSERT(ierr .eq. 0)
            call PCFactorSetFill(ctxt%pcphy, fillp, ierr)
            ASSERT(ierr .eq. 0)
            call PCFactorSetMatOrderingType(ctxt%pcphy, MATORDERINGNATURAL, ierr)
            ASSERT(ierr .eq. 0)

        else if (pre_type == mumps_pre) then
#ifdef ASTER_PETSC_HAVE_MUMPS
            ! Ou encore  mumps mais à partir du moment où petsc est compilée
            ! avec support de l'interface MUMPS
            call PCSetType(ctxt%pcphy, PCLU, ierr)
            ASSERT(ierr == 0)
            call PCFactorSetMatSolverType(ctxt%pcphy, MATSOLVERMUMPS, ierr)
            ASSERT(ierr .eq. 0)
            call PCFactorSetUpMatSolverType(ctxt%pcphy, ierr)
            ASSERT(ierr .eq. 0)
            call PCFactorGetMatrix(ctxt%pcphy, F, ierr)
            ASSERT(ierr .eq. 0)
            ! ICNTL(7) (sequential matrix ordering): 5 (METIS)
            call MatMumpsSetIcntl(F, to_petsc_int(7), to_petsc_int(5), ierr)
            ASSERT(ierr .eq. 0)
!  ICNTL(22) (in-core/out-of-core facility): 0/1
            call MatMumpsSetIcntl(F, to_petsc_int(22), to_petsc_int(1), ierr)
            ASSERT(ierr .eq. 0)
!  ICNTL(24) (detection of null pivot rows): 1
            call MatMumpsSetIcntl(F, to_petsc_int(24), to_petsc_int(1), ierr)
            ASSERT(ierr .eq. 0)
!  ICNTL(14) (percentage increase in the estimated working space)
            call MatMumpsSetIcntl(F, to_petsc_int(14), to_petsc_int(50), ierr)
            ASSERT(ierr .eq. 0)
!  CNTL(3) (absolute pivoting threshold):      1e-06
            call MatMumpsSetCntl(F, to_petsc_int(3), 1.D-6, ierr)
            ASSERT(ierr .eq. 0)
#else
            ! should not pass here
            call utmess('F', 'FERMETUR_4', sk='MUMPS')
#endif
        else
            ASSERT(.false.)
        end if
        call PCSetUp(ctxt%pcphy, ierr)
        ASSERT(ierr .eq. 0)
        !
    end subroutine set_precond_data
!
    subroutine free_augm_lagrangian_context(ctxt)
        !
        ! Dummy argument
        !
        type(augm_lagr_ctxt), intent(inout) :: ctxt
        !
        call MatDestroy(ctxt%m_mat, ierr)
        ASSERT(ierr == 0)
        call PCDestroy(ctxt%pcphy, ierr)
        ASSERT(ierr == 0)
        ! TODO compteur de référence
        call free_saddle_point_context(ctxt%sp_ctxt)
        nullify (ctxt%sp_ctxt)
        !
    end subroutine free_augm_lagrangian_context
!
#else
!
contains
!
    function new_augmented_lagrangian_context(sp_ctxt) result(ctxt)
        integer(kind=8) :: sp_ctxt
        type(augm_lagr_ctxt)         :: ctxt
        ctxt%idummy = 0
        sp_ctxt = 0
    end function new_augmented_lagrangian_context
!
    subroutine free_augm_lagrangian_context(ctxt)
        type(augm_lagr_ctxt), intent(inout) :: ctxt
        ctxt%idummy = 0
    end subroutine free_augm_lagrangian_context
#endif
end module augmented_lagrangian_context_type
