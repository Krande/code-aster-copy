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
subroutine barqua(i1, i2, coor, poin)
    implicit none
#include "asterfort/barso1.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: i1, i2, poin(*)
    real(kind=8) :: coor(*)
!     BARSOUM : TRAITEMENT DES MAILLES "QUAD8" ET "QUAD9"
!-----------------------------------------------------------------------
!
    integer(kind=8) :: i, n1, n2, n3
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!                       TRAITEMENT DES "POI1"
!     ------------------------------------------------------------------
    if (i1 .eq. 1 .and. i2 .eq. 0) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 1
                n2 = 2
                n3 = 5
            else
                n1 = 1
                n2 = 4
                n3 = 8
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
    else if (i1 .eq. 2 .and. i2 .eq. 0) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 2
                n2 = 1
                n3 = 5
            else
                n1 = 2
                n2 = 3
                n3 = 6
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
    else if (i1 .eq. 3 .and. i2 .eq. 0) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 3
                n2 = 2
                n3 = 6
            else
                n1 = 3
                n2 = 4
                n3 = 7
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
    else if (i1 .eq. 4 .and. i2 .eq. 0) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 4
                n2 = 1
                n3 = 8
            else
                n1 = 4
                n2 = 3
                n3 = 7
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
!
!     ------------------------------------------------------------------
!                       TRAITEMENT DES "SEG3"
!     ------------------------------------------------------------------
    else if (i1+i2 .eq. 3) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 1
                n2 = 4
                n3 = 8
            else
                n1 = 2
                n2 = 3
                n3 = 6
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
!
    else if (i1+i2 .eq. 7) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 4
                n2 = 1
                n3 = 8
            else
                n1 = 3
                n2 = 2
                n3 = 6
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
!
    else if ((i1+i2 .eq. 5) .and. (i1 .eq. 2 .or. i2 .eq. 2)) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 3
                n2 = 4
                n3 = 7
            else
                n1 = 2
                n2 = 1
                n3 = 5
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
!
    else if ((i1+i2 .eq. 5) .and. (i1 .eq. 4 .or. i2 .eq. 4)) then
        do i = 1, 2
            if (i .eq. 1) then
                n1 = 4
                n2 = 3
                n3 = 7
            else
                n1 = 1
                n2 = 2
                n3 = 5
            end if
            call barso1(n1, n2, n3, coor, poin)
        end do
!
    else
        call utmess('F', 'ALGORITH_36', sk='QUAD')
!
    end if
!
end subroutine
