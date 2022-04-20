! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine getHHOPara(behaviourPrepCrit)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
!
type(Behaviour_PrepCrit), intent(inout) :: behaviourPrepCrit
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get parameters from STAT_NON_LINE command
!
! --------------------------------------------------------------------------------------------------
!
! IO  behaviourPrepCrit: datastructure to prepare parameters for constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'HHO'
    integer :: hho_calc, iret, hho_stab
    character(len=16) :: txt
    real(kind=8) :: hho_coef_stab
!
! --------------------------------------------------------------------------------------------------
!
    hho_coef_stab = 0.d0
    hho_calc = HHO_CALC_NO
    hho_stab = HHO_STAB_AUTO

    call getvtx(factorKeyword, 'OPTIMISATION', iocc = 1, nbval=0, nbret = iret)

    if (iret .ne. 0) then
        if(iret == -1) then
            call getvtx(factorKeyword, 'OPTIMISATION', iocc = 1, scal=txt, nbret = iret)
            ASSERT(iret == 1)
            if(txt == "MEMOIRE") then
                hho_calc = HHO_CALC_NO
            elseif(txt == "TEMPS") then
                hho_calc = HHO_CALC_YES
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

        call getvtx(factorKeyword, 'STABILISATION', iocc = 1, nbval=0, nbret = iret)
        if(iret == -1) then
            call getvtx(factorKeyword, 'STABILISATION', iocc = 1, scal=txt, nbret = iret)
            ASSERT(iret == 1)
            if(txt == "AUTO") then
                hho_stab = HHO_STAB_AUTO
            elseif(txt == "MANUEL") then
                hho_stab = HHO_STAB_MANU
                call getvr8(factorKeyword, 'COEF_STAB', iocc = 1, scal=hho_coef_stab, nbret = iret)
                ASSERT(iret == 1)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
    end if
!
    behaviourPrepCrit%hho_coef_stab = hho_coef_stab
    behaviourPrepCrit%hho_type_stab = real(hho_stab, kind=8)
    behaviourPrepCrit%hho_type_calc = real(hho_calc, kind=8)
!
end subroutine
