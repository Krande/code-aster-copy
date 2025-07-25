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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine dbr_calcpod_savel(base, nbMode, nbSnapRedu, baseSing, baseValeR)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/romBaseSave.h"
!
    type(ROM_DS_Empi), intent(in) :: base
    integer(kind=8), intent(in) :: nbMode, nbSnapRedu
    real(kind=8), pointer :: baseValeR(:), baseSing(:)
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE - Compute
!
! Save base for lineic model
!
! --------------------------------------------------------------------------------------------------
!
! In  base             : datastructure for base
! In  nbMode           : number of modes in base
! In  nbSnapRedu       : number of snapshots used to construct base
! Ptr baseValeR        : pointer to the values of all modes in base
! Ptr baseSing         : pointer to the singular values of all modes in base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iSlice, iMode, nodeNume, cmpNume, i_2d, iEqua
    integer(kind=8) :: nbSlice, n_2d, nbCmp, nbLineMode, nbEqua
    real(kind=8), pointer :: lineModesVale(:) => null()
    real(kind=8), pointer :: lineModesSing(:) => null()
    integer(kind=8), pointer :: numeSlice(:) => null()
    type(ROM_DS_LineicNumb) :: lineicNume
!
! --------------------------------------------------------------------------------------------------
!
    lineicNume = base%lineicNume
    nbEqua = base%mode%nbEqua
    nbCmp = lineicNume%nbCmp
    nbSlice = lineicNume%nbSlice
    nbLineMode = nbMode*nbSlice
!
! - Create working objects
!
    AS_ALLOCATE(vr=lineModesVale, size=nbEqua*nbMode*nbSlice)
    AS_ALLOCATE(vr=lineModesSing, size=nbMode*nbSlice)
    AS_ALLOCATE(vi=numeSlice, size=nbMode*nbSlice)
!
! - Create index of slices
!
    do iSlice = 1, nbSlice
        do iMode = 1, nbMode
            lineModesSing(iMode+nbMode*(iSlice-1)) = baseSing(iMode)
            numeSlice(iMode+nbMode*(iSlice-1)) = iSlice
        end do
    end do
!
! - Create modes to save
!
    do iEqua = 1, nbEqua
        nodeNume = (iEqua-1)/nbCmp+1
        cmpNume = iEqua-(nodeNume-1)*nbCmp
        iSlice = lineicNume%numeSlice(nodeNume)
        n_2d = lineicNume%numeSection(nodeNume)
        i_2d = (n_2d-1)*nbCmp+cmpNume
        do iMode = 1, nbMode
            lineModesVale(iEqua+nbEqua*(iMode-1+nbMode*(iSlice-1))) = &
                baseValeR(i_2d+nbEqua/nbSlice*(iMode-1))
        end do
    end do
!
! - Save modes
!
    call romBaseSave(base, nbLineMode, nbSnapRedu, &
                     lineModesVale, lineModesSing, numeSlice)
!
! - Cleaning
!
    AS_DEALLOCATE(vr=lineModesVale)
    AS_DEALLOCATE(vr=lineModesSing)
    AS_DEALLOCATE(vi=numeSlice)
!
end subroutine
