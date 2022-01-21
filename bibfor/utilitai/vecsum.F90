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

subroutine vecsum(livr, sz, mode, vr)
  implicit none

#include "asterfort/assert.h"
#include "asterfort/twosum.h"

  real(kind=8), intent(in) :: livr(*)
  integer, intent(in) :: sz, mode
  real(kind=8), intent(out) :: vr

! ------------------------------------------------------------------------------
!
! Return the sum of the items of a vector using differents algorithms
!
! In  livr : vector of items to sum
! In  sz   : number of items to sum
! In  mode : 0 for standard summation
!            1 for Ogita-Oishi-Rump summation
!
! Out vr   : the sum value
! ------------------------------------------------------------------------------
!
  integer i
  real(kind=8) :: si, x, y

  vr = 0.0d0
  si = 0.0d0
  x  = 0.0d0
  y  = 0.0d0

  if (mode .eq. 0) then
     do i = 1, sz
        vr = vr + livr(i)
     enddo

  elseif (mode .eq. 1) then

     vr = livr(1)
     si = 0.d0
     do i = 2, sz
        call twosum(vr, livr(i), x, y)
        vr = x
        si = si + y
     enddo
     vr = vr + si

  else
     ASSERT(.false.)
  endif

end subroutine
