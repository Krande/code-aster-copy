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
subroutine maveElemCreate(base, mave_elemz, modelz)
!
    implicit none
!
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "asterfort/detrsd.h"
!
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: mave_elemz
    character(len=*), intent(in) :: modelz
!
! --------------------------------------------------------------------------------------------------
!
! Elementary vectors
!
! Create a new elementary vector
!
! --------------------------------------------------------------------------------------------------
!
! In  base             : JEVEUX base to create object
! In  mave_elem        : name of matr_elem or vect_elem
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: mave_elem, model
!
! --------------------------------------------------------------------------------------------------
!
    mave_elem = mave_elemz
    model = modelz
    call detrsd("VECT_ELEM", mave_elem)
    call vemare(base, mave_elem, model)
    call reajre(mave_elem, ' ', base)
!
end subroutine
