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
!
#include "asterf_types.h"
!
interface
    subroutine nmdcco(sddisc, i_event_acti, typdec, nbrpas, deltac,&
                      ratio , optdec      , retdec, ldcext, subdur)
        character(len=19) :: sddisc
        integer(kind=8) :: i_event_acti
        character(len=4) :: typdec
        integer(kind=8) :: nbrpas
        real(kind=8) :: deltac
        real(kind=8) :: ratio
        character(len=16) :: optdec
        integer(kind=8) :: retdec
        aster_logical :: ldcext
        real(kind=8) :: subdur
    end subroutine nmdcco
end interface
