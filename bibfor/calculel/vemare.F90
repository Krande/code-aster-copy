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

subroutine vemare(base, vect_elemz, modelz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedetr.h"
#include "asterfort/wkvect.h"
!
!
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: vect_elemz
    character(len=*), intent(in) :: modelz
!
! --------------------------------------------------------------------------------------------------
!
! RESU_ELEM Management
!
! Create RERR object for vect_elem
!
! --------------------------------------------------------------------------------------------------
!
! In  base           : JEVEUX basis
! In  matr_vect_elem : name of vect_elem
! In  modelz         : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model
    character(len=19) :: vect_elem
    character(len=24), pointer :: p_rerr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    vect_elem = vect_elemz
    model = modelz
    ASSERT(model .ne. ' ')
!
    call jedetr(vect_elem//'.RERR')
    call wkvect(vect_elem//'.RERR', base//' V K24', 1, vk24=p_rerr)
    p_rerr(1) = model
!
end subroutine
