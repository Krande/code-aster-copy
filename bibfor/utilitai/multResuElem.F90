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
subroutine multResuElem(resu_elem, coef_mult)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jaexin.h"
#include "asterfort/jexnum.h"
!
!
    character(len=19), intent(in) :: resu_elem
    real(kind=8), intent(in) :: coef_mult
!
! --------------------------------------------------------------------------------------------------
!
! Generic routines
!
! Multiply RESU_ELEM values by a constant coefficient
!
! --------------------------------------------------------------------------------------------------
!
!   In resu_elem    : name of RESU_ELEM
!   In coef_mult    : constant coefficient
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, nb_gr, nb_vale, igr
    real(kind=8), pointer :: v_vale(:) => null()
    integer(kind=8), pointer :: v_desc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

!
    call jeexin(resu_elem(1:19)//'.DESC', iret)
    ASSERT(iret .ne. 0)
    call jeveuo(resu_elem(1:19)//'.DESC', 'L', vi=v_desc)
    nb_gr = v_desc(2)

    do igr = 1, nb_gr
        call jaexin(jexnum(resu_elem(1:19)//'.RESL', igr), iret)
        if (iret .eq. 0) cycle
        call jeveuo(jexnum(resu_elem(1:19)//'.RESL', igr), 'E', vr=v_vale)
        call jelira(jexnum(resu_elem(1:19)//'.RESL', igr), 'LONMAX', nb_vale)
        v_vale(1:nb_vale) = coef_mult*v_vale(1:nb_vale)
    end do
!
end subroutine
