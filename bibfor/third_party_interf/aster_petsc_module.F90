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
! aslint: disable=C1002
module aster_petsc_module
!
    use, intrinsic :: iso_c_binding
!
#include "asterf_petsc.h"
#ifdef ASTER_HAVE_PETSC
    use petsc
!
    implicit none
!
! System routines
!
    interface
        subroutine PetscViewerAndFormatCreate(viewer, format, vf, ierr)
            use petscsysdef
            PetscViewer:: viewer
            PetscViewerFormat :: format
            PetscViewerAndFormat :: vf
            PetscErrorCode, intent(out) :: ierr
        end subroutine PetscViewerAndFormatCreate
    end interface
    interface
        subroutine PetscOptionsGetString(opt, char1, char2, char3, flg, ierr)
            use petscsysdef
            PetscOptions:: opt
            character(*) :: char1, char2, char3
            PetscBool ::  flg
            PetscErrorCode, intent(out) :: ierr
        end subroutine PetscOptionsGetString
    end interface
!
! Mat routines
!
    interface
        subroutine PCHPDDMSetAuxiliaryMat(pc, is, mat, func, ctx, ierr)
            use petsckspdef
            PC :: pc
            Mat :: mat
            IS :: is
            external :: func
            PetscInt :: ctx
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCHPDDMSetAuxiliaryMat
    end interface
    interface
        subroutine MatCreateShell(comm, m, n, mg, ng, ctxt, a_mat, ierr)
            use petscmatdef
            PetscMPIInt, intent(in) :: comm
            PetscInt, intent(in) :: m
            PetscInt, intent(in) :: n
            PetscInt, intent(in) :: mg
            PetscInt, intent(in) :: ng
            PetscInt :: ctxt
            Mat, intent(out) :: a_mat
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatCreateShell
    end interface
    ! interface
    !     subroutine MatShellSetOperation(mat, operation, myop, ierr)
    !         use petscmatdef
    !         Mat :: mat
    !         MatOperation :: operation
    !         external :: myop
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine MatShellSetOperation
    ! end interface
!
! PC and KSP routines
!
    ! interface
    !     subroutine PCFactorSetMatOrderingType(pc, ordering, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         character(*) :: ordering
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCFactorSetMatOrderingType
    ! end interface
    ! interface
    !     subroutine PCShellSetSetup(pc, mysetup, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         external :: mysetup
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetSetUp
    ! end interface
    ! interface
    !     subroutine PCShellSetApply(pc, myapply, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         external :: myapply
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetApply
    ! end interface
    ! interface
    !     subroutine PCShellSetApplySymmetricRight(pc, myapply, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         external :: myapply
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetApplySymmetricRight
    ! end interface
    ! interface
    !     subroutine PCShellSetApplySymmetricLeft(pc, myapply, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         external :: myapply
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetApplySymmetricLeft
    ! end interface
    ! interface
    !     subroutine PCShellSetDestroy(pc, mydestroy, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         external :: mydestroy
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetDestroy
    ! end interface
    interface
        subroutine KSPMonitorSet(ksp, mykspmonitor, vf, mydestroy, ierr)
            use petsckspdef
            KSP :: ksp
            PetscViewerAndFormat:: vf
            external :: mykspmonitor, mydestroy
            PetscErrorCode, intent(out) :: ierr
        end subroutine KSPMonitorSet
    end interface
    ! interface
    !     subroutine PCShellSetName(pc, myname, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         character(*) :: myname
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCShellSetName
    ! end interface
    ! interface
    !     subroutine PCFactorSetMatSolverType(pc, solvertype, ierr)
    !         use petsckspdef
    !         PC :: pc
    !         character(*):: solvertype
    !         PetscErrorCode, intent(out) :: ierr
    !     end subroutine PCFactorSetMatSolverType
    ! end interface
#endif
end module aster_petsc_module
