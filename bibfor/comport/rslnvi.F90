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

subroutine rslnvi(elem_model, ndt_, ndi_, nr_, nvi_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
!
    character(len=*), intent(in) :: elem_model
    integer(kind=8), optional, intent(out) :: ndt_
    integer(kind=8), optional, intent(out) :: ndi_
    integer(kind=8), optional, intent(out) :: nr_
    integer(kind=8), optional, intent(out) :: nvi_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment ROUSSELIER
!
! Size of tensors and number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_model : type of model on element
! Out ndt        : number of terms in tensors
! Out ndi        : number of "indirect" (?) terms in tensors
! Out nr         : number of non-linear system
! Out nvi        : number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndt, ndi, nr, nvi
!
! --------------------------------------------------------------------------------------------------
!
    nvi = 5
    if (elem_model(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
        nr = ndt+2
    else if (elem_model(1:6) .eq. 'D_PLAN' .or. elem_model(1:4) .eq. 'AXIS') then
        ndt = 4
        ndi = 3
        nr = ndt+2
    else if (elem_model(1:6) .eq. 'C_PLAN') then
        call utmess('F', 'COMPOR5_52')
    else if (elem_model(1:2) .eq. '1D') then
        ndt = 3
        ndi = 3
        nr = ndt+2
    else
        ASSERT(.false.)
    end if
!
    if (present(ndt_)) then
        ndt_ = ndt
    end if
    if (present(ndi_)) then
        ndi_ = ndi
    end if
    if (present(nr_)) then
        nr_ = nr
    end if
    if (present(nvi_)) then
        nvi_ = nvi
    end if
!
end subroutine
