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
!
subroutine te0125(option, nomte)
!
    use SolidShell_Elementary_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_SOLIDE (HEXA9 and PENTA7)
!
! Options: RIGI_MECA, MASS_MECA, RIGI_GEOM
!          SIEF_ELGA, FORC_NODA, EPSI_ELGA, EPSL_ELGA
!          CHAR_MECA_PRES_R, CHAR_MECA_FF3D3D, CHAR_MECA_FR3D3D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    if (option .eq. 'RIGI_MECA') then
        call compRigiMatr()
    elseif (option .eq. 'SIEF_ELGA') then
        call compSiefElga()
    elseif (option .eq. 'FORC_NODA') then
        call compForcNoda()
    elseif (option .eq. 'REFE_FORC_NODA') then
        call compRefeForcNoda()
    elseif (option .eq. 'EPSI_ELGA') then
        call compEpsiElga()
    elseif (option .eq. 'EPSL_ELGA') then
        call compEpslElga()
    elseif (option .eq. 'EPVC_ELGA') then
        call compEpvcElga()
    elseif (option .eq. 'CHAR_MECA_PRES_R') then
        call compLoad(option)
    elseif (option .eq. 'CHAR_MECA_PESA_R') then
        call compLoad(option)
    elseif (option .eq. 'CHAR_MECA_FF3D3D') then
        call compLoad(option)
    elseif (option .eq. 'CHAR_MECA_FR3D3D') then
        call compLoad(option)
    elseif (option .eq. 'CHAR_MECA_TEMP_R') then
        call compLoadExteStatVari(option)
    elseif (option .eq. 'CHAR_MECA_HYDR_R') then
        call compLoadExteStatVari(option)
    elseif (option .eq. 'CHAR_MECA_SECH_R') then
        call compLoadExteStatVari(option)
    elseif (option .eq. 'CHAR_MECA_EPSA_R') then
        call compLoadExteStatVari(option)
    elseif (option .eq. 'MASS_MECA') then
        call compMassMatr()
    elseif (option .eq. 'RIGI_GEOM') then
        call compRigiGeomMatr()
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
