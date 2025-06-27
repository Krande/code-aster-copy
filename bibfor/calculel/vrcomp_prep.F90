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
subroutine vrcomp_prep(variZ, variRedu, &
                       comporCurrZ, comporCurrRedu, &
                       comporPrevZ_, comporPrevRedu_)
!
    implicit none
!
#include "asterfort/carces.h"
#include "asterfort/cesred.h"
#include "asterfort/cestas.h"
#include "asterfort/celces.h"
#include "asterfort/detrsd.h"
!
!
    character(len=*), intent(in) :: variZ
    character(len=19), intent(out) :: variRedu
    character(len=*), intent(in)  :: comporCurrZ
    character(len=19), intent(out) :: comporCurrRedu
    character(len=*), optional, intent(in)  :: comporPrevZ_
    character(len=19), optional, intent(out) :: comporPrevRedu_
!
! --------------------------------------------------------------------------------------------------
!
! Check compatibility of comportments
!
! Prepare fields
!
! --------------------------------------------------------------------------------------------------
!
! In  vari           : internal variable
! Out variRedu       : reduced field for internal variable
! In  comporCurr     : current comportment
! Out comporCurrRedu : reduced field for current comportment
! In  comporPrev     : previous comportment
! Out comporPrevRedu : reduced field for previous comportment
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: coto = '&&VRCOMP.COTO'
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!

! - Create reduced CARTE on current comportement
    comporCurrRedu = '&&VRCOMP.COPP'
    call carces(comporCurrZ, 'ELEM', ' ', 'V', coto, 'A', iret)
    call cesred(coto, 0, [0], 1, 'RELCOM', 'V', comporCurrRedu)
    call detrsd('CHAM_ELEM_S', coto)

! - Create reduced field for internal state variables
    variRedu = '&&VRCOMP.VARI_R'
    call celces(variZ, 'V', variRedu)
    call cestas(variRedu)

! - Create reduced CARTE on previous comportement
    if (present(comporPrevZ_)) then
        comporPrevRedu_ = '&&VRCOMP.COPM'
        call carces(comporPrevZ_, 'ELEM', ' ', 'V', coto, 'A', iret)
        call cesred(coto, 0, [0], 1, 'RELCOM', 'V', comporPrevRedu_)
    end if

    call detrsd('CHAM_ELEM_S', coto)
!
end subroutine
