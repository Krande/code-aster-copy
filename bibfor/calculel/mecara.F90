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
subroutine mecara(caraElem, chcara)
!
    implicit none
!
#include "asterf_types.h"
!
    character(len=*), intent(in) :: caraElem
    character(len=*), intent(inout) :: chcara(18)
!
! --------------------------------------------------------------------------------------------------
!
! Preparing <CARTE> field for elementary characteristics
!
! --------------------------------------------------------------------------------------------------
!
! In  caraElem    : name of elementary characteristics (field)
! IO  chcara      : name of <CARTE> field for elementary characteristics
!
! --------------------------------------------------------------------------------------------------
!
    chcara = ' '
!
    if (caraElem(1:8) .ne. ' ') then
        chcara(1) = caraElem(1:8)//'.CARORIEN'
        chcara(2) = caraElem(1:8)//'.CARDISCK'
        chcara(3) = caraElem(1:8)//'.CARDISCM'
        chcara(4) = caraElem(1:8)//'.CARDISCA'
        chcara(5) = caraElem(1:8)//'.CARGEOPO'
        chcara(6) = caraElem(1:8)//'.CARGENPO'
        chcara(7) = caraElem(1:8)//'.CARCOQUE'
        chcara(9) = caraElem(1:8)//'.CARARCPO'
        chcara(10) = caraElem(1:8)//'.CARCABLE'
        chcara(11) = caraElem(1:8)//'.CARGENBA'
        chcara(12) = caraElem(1:8)//'.CARMASSI'
        chcara(13) = caraElem(1:8)//'.CARPOUFL'
        chcara(14) = caraElem(1:8)//'.CVENTCXF'
        chcara(15) = caraElem(1:8)//'.CARDINFO'
        chcara(16) = caraElem(1:8)//'.CANBSP'
        chcara(17) = caraElem(1:8)//'.CAFIBR'
    end if
!
end subroutine
