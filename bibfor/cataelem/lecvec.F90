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
subroutine lecvec(iad, long, type, unite)
    implicit none
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
    character(len=3) :: type
    integer(kind=8) :: unite
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iad, k, long
!-----------------------------------------------------------------------
    if (type .eq. 'R') then
        do k = 1, long
            read (unite, '(1E12.5)') zr(iad-1+k)
        end do
!
    else if (type .eq. 'I') then
        do k = 1, long
            read (unite, '(I12)') zi(iad-1+k)
        end do
!
    else if (type .eq. 'K8') then
        do k = 1, long
            read (unite, '(A8)') zk8(iad-1+k)
        end do
!
    else if (type .eq. 'K16') then
        do k = 1, long
            read (unite, '(A16)') zk16(iad-1+k)
        end do
!
    else if (type .eq. 'K24') then
        do k = 1, long
            read (unite, '(A24)') zk24(iad-1+k)
        end do
!
    else if (type .eq. 'K32') then
        do k = 1, long
            read (unite, '(A32)') zk32(iad-1+k)
        end do
!
    else if (type .eq. 'K80') then
        do k = 1, long
            read (unite, '(A80)') zk80(iad-1+k)
        end do
    else
        ASSERT(.false.)
    end if
!
!
end subroutine
