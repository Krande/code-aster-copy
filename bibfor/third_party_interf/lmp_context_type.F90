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
! A lmp_context object is a container used to manage
! a Limited Memory Preconditioner (second level preconditioner)
!
! It contains
! - parameters of the method : ritz
! - the necessary data to apply the preconditioner to a PETSc Vec array
!
module lmp_context_type
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
! person_in_charge: natacha.bereux at edf.fr
!
    use aster_petsc_module
    implicit none
!
    private
#include "asterf.h"
#include "asterfort/assert.h"
!
!
    type, public :: lmp_ctxt
        !
#ifdef ASTER_HAVE_PETSC
        !
        ! Nombre de vecteurs de Ritz souhaité pour la construction du préconditionneur
        integer(kind=8)  :: ritz
        ! Nombre de vecteurs de Ritz réellement utilisés dans chaque LMP (vaut ritz ou ritz+1)
        PetscInt :: ritzeff
        ! ss contient les vecteurs propres de Ritz
        Vec, dimension(:), allocatable   :: ss
        ! yy et zz servent à appliquer le LMP
        Vec, dimension(:), allocatable  :: yy, zz
        ! Valeurs propres de Ritz  (parties réelles et imaginaires)
        PetscReal, dimension(:), allocatable :: thetar, thetai
        ! Vecteur de travail
        PetscReal, dimension(:), allocatable :: sk
#endif
    end type lmp_ctxt
!
    public :: new_lmp_context, build_lmp_context, free_lmp_context
!
#ifdef ASTER_HAVE_PETSC
    PetscErrorCode :: ierr
!
contains
!
! Allocate a fresh new context
!
    function new_lmp_context(ritz) result(ctxt)
        !
        integer(kind=8), intent(in), optional   :: ritz
        type(lmp_ctxt)          :: ctxt
        !
        ! Local variables
        !
        integer(kind=8) :: ierr
        !
        if (present(ritz)) then
            ctxt%ritz = ritz
        else
            ctxt%ritz = 5
        end if
        !
        ctxt%ritzeff = to_petsc_int(ctxt%ritz)
        !
        allocate (ctxt%ss(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        allocate (ctxt%yy(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        allocate (ctxt%zz(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        allocate (ctxt%thetar(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        allocate (ctxt%thetai(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        allocate (ctxt%sk(ctxt%ritz+1), stat=ierr)
        ASSERT(ierr == 0)
        !
    end function new_lmp_context
!
!
! Init data needed by a LMP  preconditioner
! On entry KSP object
! It is mandatory to :
! - build the lmp context ONLY after a first solve with the ksp
! - to setup the Ritz computation in the ksp before this solve
!
    subroutine build_lmp_context(ksp, ctxt)
        !
        ! Dummy arguments
        !
        KSP, intent(in) :: ksp
        type(lmp_ctxt), intent(inout) :: ctxt
        !
        ! Local variables
        !
        PetscErrorCode :: ierr
        Vec :: xx, sn(3)
        PetscInt :: ii, jj, mm
        PC :: pc
        Mat :: amat
        PetscScalar :: norm
        !
        call KSPGetOperators(ksp, amat, PETSC_NULL_MAT, ierr)
        ASSERT(ierr == 0)
        ! On calcule les vecteurs de Ritz de plus petit module pour definir le LMP
        ! allocation des vecteurs PETSc à partir du vecteur solution du ksp
        call KSPGetSolution(ksp, xx, ierr)
        ASSERT(ierr == 0)
        mm = ctxt%ritz+1
        call VecDuplicateVecs(xx, mm, ctxt%ss, ierr)
        ASSERT(ierr == 0)
        ! Calcul de ritzeff vecteurs de Ritz
        call KSPComputeRitz(ksp, PETSC_TRUE, PETSC_TRUE, ctxt%ritzeff, ctxt%ss, &
                            ctxt%thetar, ctxt%thetai, ierr)
        ASSERT(ierr == 0)
        !
!   On doit encore effectuer une orthonormalisation pour faciliter l'application du
!   préconditionneur
!       H=I-KS(S^TK^TKS)S^TK^T+S(S^TK^TKS)S^TK^T
!   =>  H=I+ZY^T
!
!   Creation des vecteurs de travail sn
        mm = 3
        call VecDuplicateVecs(ctxt%ss(1), mm, sn, ierr)
        ASSERT(ierr == 0)

        call KSPGetPC(ksp, pc, ierr)
        ASSERT(ierr == 0)
!
!   Creation des matrices qui stockent les directions A^{T}A-orthonormales
        mm = ctxt%ritz+1
        call VecDuplicateVecs(ctxt%ss(1), mm, ctxt%yy, ierr)
        ASSERT(ierr == 0)
        call VecDuplicateVecs(ctxt%ss(1), mm, ctxt%zz, ierr)
        ASSERT(ierr == 0)

!   Algorithme de Gram Schmidt (thèse de Sylvain Mercier p. 111, algo 12)
!
!   Iteration 0
        call VecCopy(ctxt%ss(1), sn(1), ierr)
        ASSERT(ierr == 0)
        call PCApplySymmetricRight(pc, sn(1), sn(2), ierr)
        ASSERT(ierr == 0)
        call MatMult(amat, sn(2), sn(3), ierr)
        ASSERT(ierr == 0)
        call PCApplySymmetricLeft(pc, sn(3), sn(2), ierr)
        ASSERT(ierr == 0)
        call VecNorm(sn(2), NORM_2, norm, ierr)
        ASSERT(ierr == 0)
        call VecScale(sn(1), 1.d0/norm, ierr)
        ASSERT(ierr == 0)
        call VecScale(sn(2), 1.d0/norm, ierr)
        ASSERT(ierr == 0)
        call VecCopy(sn(1), ctxt%zz(1), ierr)
        ASSERT(ierr == 0)
        call VecCopy(sn(2), ctxt%yy(1), ierr)
        ASSERT(ierr == 0)
!
!   Boucle sur le nombre de directions
        do ii = 1, ctxt%ritzeff-1
            call VecCopy(ctxt%ss(ii+1), sn(1), ierr)
            ASSERT(ierr == 0)
            call PCApplySymmetricRight(pc, sn(1), sn(2), ierr)
            ASSERT(ierr == 0)
            call MatMult(amat, sn(2), sn(3), ierr)
            ASSERT(ierr == 0)
            call PCApplySymmetricLeft(pc, sn(3), sn(2), ierr)
            ASSERT(ierr == 0)
            do jj = 1, ii
                call VecDot(ctxt%yy(jj), sn(2), ctxt%sk(jj), ierr)
                ASSERT(ierr == 0)
                ctxt%sk(jj) = -1.d0*ctxt%sk(jj)
            end do
            call VecMAXPY(sn(1), ii, ctxt%sk, ctxt%zz, ierr)
            ASSERT(ierr == 0)
            call VecMAXPY(sn(2), ii, ctxt%sk, ctxt%yy, ierr)
            ASSERT(ierr == 0)
            call VecNorm(sn(2), NORM_2, norm, ierr)
            ASSERT(ierr == 0)
            call VecScale(sn(1), 1.d0/norm, ierr)
            ASSERT(ierr == 0)
            call VecScale(sn(2), 1.d0/norm, ierr)
            ASSERT(ierr == 0)
            call VecCopy(sn(1), ctxt%zz(ii+1), ierr)
            ASSERT(ierr == 0)
            call VecCopy(sn(2), ctxt%yy(ii+1), ierr)
            ASSERT(ierr == 0)
        end do
        do ii = 1, ctxt%ritzeff
            call VecAXPY(ctxt%zz(ii), -1.d0, ctxt%yy(ii), ierr)
            ASSERT(ierr == 0)
        end do
        !
!   Libérations mémoire
        do ii = 1, 3
            call VecDestroy(sn(ii), ierr)
            ASSERT(ierr == 0)
        end do
!
    end subroutine build_lmp_context
!
!
    subroutine free_lmp_context(ctxt)
        !
        ! Dummy argument
        !
        type(lmp_ctxt), intent(inout) :: ctxt
        ! Local Variables
        integer(kind=8) :: ii
        !
        if (allocated(ctxt%ss)) then
            do ii = 1, size(ctxt%ss)
                call VecDestroy(ctxt%ss(ii), ierr)
                ASSERT(ierr == 0)
            end do
            deallocate (ctxt%ss)
        end if
        if (allocated(ctxt%yy)) then
            do ii = 1, size(ctxt%yy)
                call VecDestroy(ctxt%yy(ii), ierr)
                ASSERT(ierr == 0)
            end do
            deallocate (ctxt%yy)
        end if
        if (allocated(ctxt%zz)) then
            do ii = 1, size(ctxt%zz)
                call VecDestroy(ctxt%zz(ii), ierr)
                ASSERT(ierr == 0)
            end do
            deallocate (ctxt%zz)
        end if
        if (allocated(ctxt%thetar)) then
            deallocate (ctxt%thetar)
        end if
        if (allocated(ctxt%thetai)) then
            deallocate (ctxt%thetai)
        end if
        if (allocated(ctxt%sk)) then
            deallocate (ctxt%sk)
        end if
        !
    end subroutine free_lmp_context
!
#else
!
contains
!
    function new_lmp_context(ritz) result(ctxt)
        integer(kind=8)  :: ritz
        integer(kind=8)  :: ctxt
        ritz = 0
        ctxt = 0
        ASSERT(.false.)
    end function new_lmp_context
!
    subroutine build_lmp_context(ksp, ctxt)
        integer(kind=8) :: ksp
        integer(kind=8) :: ctxt
        ksp = 0
        ctxt = 0
        ASSERT(.false.)
    end subroutine build_lmp_context
!
    subroutine free_lmp_context(ctxt)
        integer(kind=8) :: ctxt
        ctxt = 0
        ASSERT(.false.)
    end subroutine free_lmp_context
#endif
end module lmp_context_type
