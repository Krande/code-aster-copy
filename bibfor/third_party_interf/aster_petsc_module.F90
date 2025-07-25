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
#include "asterf_petsc.h"
#ifdef ASTER_HAVE_PETSC
    use petscsysdef
    use petscvecdef
    use petscmatdef
    use petscpcdef
    use petsckspdef
    use petscsys
    use petscvec
    use petscmat
    use petscpc
    use petscksp
!
    implicit none
!
    interface
        subroutine PetscObjectSetName(obj, description, ierr)
            use petscsysdef
            type(*) :: obj
            character(len=*) :: description
            PetscErrorCode, intent(out) :: ierr
        end subroutine PetscObjectSetName
    end interface
    interface
        subroutine PetscLogDefaultBegin(ierr)
            use petscsysdef
            PetscErrorCode, intent(out) :: ierr
        end subroutine PetscLogDefaultBegin
    end interface
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
        subroutine PetscObjectGetComm(vec, comm, ierr)
            use petscvecdef
            Vec, intent(in):: vec
            PetscMPIInt, intent(out) :: comm
            PetscErrorCode, intent(out) :: ierr
        end subroutine PetscObjectGetComm
    end interface
!
! Vec routines
!
    interface
        subroutine VecDuplicateVecs(v, m, vtab, ierr)
            use petscvecdef
            Vec, intent(in)  :: v
            PetscInt, intent(in) :: m
            Vec  :: vtab(*)
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecDuplicateVecs
    end interface
    interface
        subroutine VecGetArray(x, x_array, i_x, ierr)
            use petscvecdef
            Vec :: x
            PetscScalar :: x_array(*)
            PetscOffset :: i_x
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecGetArray
    end interface
    interface
        subroutine VecRestoreArray(x, x_array, i_x, ierr)
            use petscvecdef
            Vec :: x
            PetscScalar :: x_array(*)
            PetscOffset :: i_x
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecRestoreArray
    end interface
!
! VecScatter routines
!
    interface
        subroutine VecScatterCreate(x, ix, y, iy, newsf, ierr)
            use petscvecdef
            Vec, intent(in) :: x, y
            IS, intent(in) :: ix, iy
            VecScatter, intent(out)  :: newsf
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecScatterCreate
    end interface
    interface
        subroutine VecScatterDestroy(v, ierr)
            use petscvecdef
            VecScatter, intent(in)  :: v
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecScatterDestroy
    end interface
    interface
        subroutine VecScatterBegin(sf, x, y, addv, mode, ierr)
            use petscvecdef
            Vec, intent(in) :: x, y
            InsertMode, intent(in) :: addv
            ScatterMode, intent(in) :: mode
            VecScatter, intent(in)  :: sf
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecScatterBegin
    end interface
    interface
        subroutine VecScatterEnd(sf, x, y, addv, mode, ierr)
            use petscvecdef
            Vec, intent(in) :: x, y
            InsertMode, intent(in) :: addv
            ScatterMode, intent(in) :: mode
            VecScatter, intent(in)  :: sf
            PetscErrorCode, intent(out) :: ierr
        end subroutine VecScatterEnd
    end interface
!
! Mat routines
!
    interface
        subroutine MatConvert(mat, newtype, reuse, matnew, ierr)
            use petscmatdef
            Mat :: mat
            character(*) :: newtype
            MatReuse :: reuse
            Mat :: matnew
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatConvert
    end interface
    interface
        subroutine PCHPDDMSetAuxiliaryMat(pc, is, mat, func, ctx, ierr)
            use petsckspdef
            PC :: pc
            Mat :: mat
            IS :: is
            external :: func
            PetscInt :: ctx(*)
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
            PetscInt, intent(in) :: ctxt(*)
            Mat, intent(out) :: a_mat
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatCreateShell
    end interface
    interface
        subroutine MatCreateSubMatrices(mat, n, irow, icol, scall, submat, ierr)
            use petscmatdef
            Mat :: mat
            PetscInt :: n
            IS :: irow, icol
            MatReuse :: scall
            Mat :: submat(*)
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatCreateSubMatrices
    end interface
    interface
        subroutine MatCreateVecs(a, vright, vleft, ierr)
            use petscmatdef
            Mat, intent(in) :: a
            Vec :: vright, vleft
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatCreateVecs
    end interface
    interface
        subroutine MatGetRowIJ(mat, shift, symmetric, inodecompressed, n, ia, iia, &
                               ja, jja, done, ierr)
            use petscmatdef
            Mat :: mat
            PetscInt :: shift
            PetscBool :: symmetric, inodecompressed
            PetscInt :: n
            PetscInt :: ia(*), ja(*)
            PetscOffset :: iia, jja
            PetscBool :: done
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatGetRowIJ
    end interface
    interface
        subroutine MatRestoreRowIJ(mat, shift, symmetric, inodecompressed, n, ia, iia, &
                                   ja, jja, done, ierr)
            use petscmatdef
            Mat :: mat
            PetscInt :: shift
            PetscBool :: symmetric, inodecompressed
            PetscInt :: n
            PetscInt :: ia(*), ja(*)
            PetscOffset :: iia, jja
            PetscBool :: done
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatRestoreRowIJ
    end interface
    interface
        subroutine MatShellSetOperation(mat, operation, myop, ierr)
            use petscmatdef
            Mat :: mat
            MatOperation :: operation
            external :: myop
            PetscErrorCode, intent(out) :: ierr
        end subroutine MatShellSetOperation
    end interface
!
! PC and KSP routines
!
    interface
        subroutine PCFactorSetMatOrderingType(pc, ordering, ierr)
            use petsckspdef
            PC :: pc
            character(*) :: ordering
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCFactorSetMatOrderingType
    end interface
    interface
        subroutine PCFieldSplitGetSubKSP(pc, nsplit, subksp, ierr)
            use petsckspdef
            PC :: pc
            PetscInt :: nsplit
            KSP  :: subksp(*)
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCFieldSplitGetSubKSP
    end interface
    interface
        subroutine PCShellSetSetup(pc, mysetup, ierr)
            use petsckspdef
            PC :: pc
            external :: mysetup
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetSetUp
    end interface
    interface
        subroutine PCShellSetApply(pc, myapply, ierr)
            use petsckspdef
            PC :: pc
            external :: myapply
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetApply
    end interface
    interface
        subroutine PCShellSetApplySymmetricRight(pc, myapply, ierr)
            use petsckspdef
            PC :: pc
            external :: myapply
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetApplySymmetricRight
    end interface
    interface
        subroutine PCShellSetApplySymmetricLeft(pc, myapply, ierr)
            use petsckspdef
            PC :: pc
            external :: myapply
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetApplySymmetricLeft
    end interface
    interface
        subroutine PCShellSetDestroy(pc, mydestroy, ierr)
            use petsckspdef
            PC :: pc
            external :: mydestroy
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetDestroy
    end interface
    interface
        subroutine PCShellSetName(pc, myname, ierr)
            use petsckspdef
            PC :: pc
            character(*) :: myname
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCShellSetName
    end interface
    interface
        subroutine PCFactorSetMatSolverType(pc, solvertype, ierr)
            use petsckspdef
            PC :: pc
            character(*):: solvertype
            PetscErrorCode, intent(out) :: ierr
        end subroutine PCFactorSetMatSolverType
    end interface
    interface
        subroutine KSPMonitorSet(ksp, mykspmonitor, vf, mydestroy, ierr)
            use petsckspdef
            KSP :: ksp
            PetscViewerAndFormat:: vf
            external :: mykspmonitor, mydestroy
            PetscErrorCode, intent(out) :: ierr
        end subroutine KSPMonitorSet
    end interface
#endif
end module aster_petsc_module
