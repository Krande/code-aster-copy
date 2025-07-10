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
subroutine getTHMPara(prepMapCarcri)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/getvr8.h"
!
    type(BehaviourPrep_MapCarcri), intent(inout) :: prepMapCarcri
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get parameters from STAT_NON_LINE command
!
! --------------------------------------------------------------------------------------------------
!
! IO  prepMapCarcri    : datastructure to construct CARCRI map
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'SCHEMA_THM'
    integer(kind=8) :: iret
    real(kind=8) :: parm_theta_thm, parm_alpha_thm
!
! --------------------------------------------------------------------------------------------------
!
    parm_theta_thm = 1.d0
    parm_alpha_thm = 1.d0

    call getvr8(factorKeyword, 'PARM_THETA', iocc=1, nbval=0, nbret=iret)
    if (iret .ne. 0) then
        call getvr8(factorKeyword, 'PARM_THETA', iocc=1, scal=parm_theta_thm, nbret=iret)
        call getvr8(factorKeyword, 'PARM_ALPHA', iocc=1, scal=parm_alpha_thm, nbret=iret)
    end if
!
    prepMapCarcri%parm_theta_thm = parm_theta_thm
    prepMapCarcri%parm_alpha_thm = parm_alpha_thm
!
end subroutine
