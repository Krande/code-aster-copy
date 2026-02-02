! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine metaGetType(metaType, nbPhases)
!
    use calcul_module, only: calcul_status
    implicit none
!
#include "asterfort/rcvarc.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(out) :: metaType, nbPhases
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get metallurgy type
!
! --------------------------------------------------------------------------------------------------
!
! Out metaType       : type of metallurgy
! Out nbPhases       : total number of phases (cold and hot)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpg = 1, ksp = 1
    character(len=8), parameter :: steel = "PFERRITE", zirc = "ALPHPUR"
    character(len=8) :: fami
    integer(kind=8) :: iret_steel, iret_zirc
    real(kind=8) :: r8dummy
!
! --------------------------------------------------------------------------------------------------
!
    metaType = META_NONE
    nbPhases = 0

! - Choice of integration scheme: for CALC_POINT_MAT is PMAT !
    if (calcul_status() .eq. 2) then
        fami = 'PMAT'
    else
        fami = 'RIGI'
    end if
!
    call rcvarc(' ', steel, '+', fami, kpg, &
                ksp, r8dummy, iret_steel)
    if (iret_steel .eq. 0) then
        metaType = META_STEEL
        nbPhases = 5

    else
        call rcvarc(' ', zirc, '+', fami, kpg, &
                    ksp, r8dummy, iret_zirc)
        if (iret_zirc .eq. 0) then
            metaType = META_ZIRC
            nbPhases = 3
        else
            metaType = META_NONE
            nbPhases = 0
        end if

    end if
!
end subroutine
