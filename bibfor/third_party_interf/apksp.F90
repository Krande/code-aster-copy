! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

subroutine apksp(kptsc)
!
#include "asterf_types.h"
#include "asterf_petsc.h"

    use aster_petsc_module
    use petsc_data_module
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  CREATION DU SOLVEUR DE KRYLOV PETSC (INSTANCE NUMERO KPTSC)
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: nmaxit, ifm, niv
!
    character(len=24) :: algo, precon
    character(len=19) :: nomat, nosolv
    character(len=14) :: nonu
!
    real(kind=8) :: resire
    real(kind=8), pointer :: slvr(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscErrorCode :: ierr
    PetscInt :: maxits
    PetscReal :: rtol, atol, dtol, aster_petsc_real
    Mat :: a
    KSP :: ksp
    PetscViewerAndFormat :: vf
!=================================================================
    call jemarq()
!
    call infniv(ifm, niv)
!
    aster_petsc_real = PETSC_DEFAULT_REAL
!     -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
    ksp = kp(kptsc)
!
    call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
    call jeveuo(nosolv//'.SLVR', 'L', vr=slvr)
    call jeveuo(nosolv//'.SLVI', 'L', vi=slvi)
    algo = slvk(6)
    precon = slvk(2)
    resire = slvr(2)
    nmaxit = slvi(2)
!
!     -- choix de l'algo KSP :
!     ------------------------
    if (precon .ne. 'UTILISATEUR') then
        if (algo .eq. 'CG') then
            call KSPSetType(ksp, KSPCG, ierr)
            ASSERT(ierr .eq. 0)
        else if (algo .eq. 'CR') then
            call KSPSetType(ksp, KSPCR, ierr)
            ASSERT(ierr .eq. 0)
        else if (algo .eq. 'GMRES') then
            call KSPSetType(ksp, KSPGMRES, ierr)
            call KSPGMRESSetBreakdownTolerance(ksp, 1.d3, ierr)
            ASSERT(ierr .eq. 0)
        else if (algo .eq. 'GMRES_LMP') then
            call KSPSetType(ksp, KSPGMRES, ierr)
            call KSPGMRESSetBreakdownTolerance(ksp, 1.d3, ierr)
            ASSERT(ierr .eq. 0)
            ! On indique qu'on veut récupérer les Ritz à la fin de la résolution
            ! afin de construire le préconditionneur de second niveau LMP
            call KSPSetComputeRitz(ksp, petsc_true, ierr)
            ASSERT(ierr .eq. 0)
        else if (algo .eq. 'GCR') then
            call KSPSetType(ksp, KSPGCR, ierr)
            ASSERT(ierr .eq. 0)
        else if (algo .eq. 'FGMRES') then
            call KSPSetType(ksp, KSPFGMRES, ierr)
            ASSERT(ierr .eq. 0)
        else
            ASSERT(.false.)
        end if
    end if
!
!     -- paramètres numériques :
!     --------------------------
!
!     -- nb iter max :
    if (nmaxit .le. 0) then
        maxits = PETSC_DEFAULT_INTEGER
    else
        maxits = to_petsc_int(nmaxit)
    end if
!
    rtol = resire
    atol = aster_petsc_real
    dtol = aster_petsc_real
!
    call KSPSetTolerances(ksp, rtol, atol, dtol, maxits, ierr)
    ASSERT(ierr .eq. 0)
!   on insert les options utilisateur juste avant l'appel à xxSetFromOptions
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
    ASSERT(ierr .eq. 0)
    call KSPSetFromOptions(ksp, ierr)
    ASSERT(ierr == 0)
    call KSPSetUp(ksp, ierr)
    ASSERT(ierr == 0)
!
!     - pour suivre les itérations de Krylov
!     --------------------------------------
    if (niv .ge. 2) then
        call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr)
        call KSPMonitorSet(ksp, KSPMonitorTrueResidual, vf, PetscViewerAndFormatDestroy, &
                           ierr)
        ASSERT(ierr .eq. 0)
    end if
!
    call jedema()
!
#else
    kptsc = 0
    ASSERT(.false.)
#endif
!
end subroutine
