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

subroutine memare(base, matr_vect_elemz, modelz, suropt, l_ss)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedetr.h"
#include "asterfort/wkvect.h"
!
!
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: matr_vect_elemz
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: suropt
    aster_logical, optional, intent(in) :: l_ss
!
! --------------------------------------------------------------------------------------------------
!
! RESU_ELEM Management
!
! Create RERR object for matr_elem or vect_elem
!
! --------------------------------------------------------------------------------------------------
!
! In  base           : JEVEUX basis
! In  matr_vect_elem : name of matr_elem or vect_elem
! In  modelz         : name of model
! In  suropt         : name of "SUR_OPTION"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model
    character(len=19) :: matr_vect_elem
    character(len=24), pointer :: p_rerr(:) => null()
    aster_logical :: lss
!
! --------------------------------------------------------------------------------------------------
!
    matr_vect_elem = matr_vect_elemz
    model = modelz
    ASSERT(model .ne. ' ')
!
    call jedetr(matr_vect_elem//'.RERR')
    call wkvect(matr_vect_elem//'.RERR', base//' V K24', 3, vk24=p_rerr)
    p_rerr(1) = model
    p_rerr(2) = suropt
    p_rerr(3) = "NON_SOUS_STRUC"
!
    lss = ASTER_FALSE
    if (present(l_ss)) then
        lss = l_ss
    end if
    if (lss) then
        p_rerr(3) = "OUI_SOUS_STRUC"
    end if
!
end subroutine
