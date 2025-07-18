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
subroutine num2k8(nomgd, tglok8, tlock8, nblk8, tind)
    implicit none
    integer(kind=8) :: nblk8, tind(*)
    character(len=8) :: nomgd, tglok8(*), tlock8(*)
!     COPIE DE NUMEK8
!     TIND(I) <-- INDICE DANS LE TABLEAU TGLOK8 DE L' ELEMEMT
!                 NUMERO I DE TLOCK8
!                 (NBLK8 : DIMENSION DE TLOCK8)
!
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
!-----------------------------------------------------------------------
    if (nomgd(1:6) .eq. 'SIEF_R') then
        do i = 1, nblk8, 1
            tind(i) = 0
!  COMPOSANTES TRAITEES: SIXX SIYY SIZZ SIXY SIXZ SIYZ
            do j = 1, 6, 1
                if (tlock8(i) .eq. tglok8(j)) then
                    tind(i) = j
                    goto 10
                end if
            end do
!  COMPOSANTES TRAITEES: NXX NYY NXY MXX MYY MXY
            do j = 14, 19, 1
                if (tlock8(i) .eq. tglok8(j)) then
                    tind(i) = j
                    goto 10
                end if
            end do
10          continue
        end do
!
    else if (nomgd(1:6) .eq. 'EPSI_R') then
!  COMPOSANTES TRAITEES: EPXX EPYY EPZZ EPXY EPXZ EPYZ
!                        EXX EYY EXY KXX KYY KXY
        do i = 1, nblk8, 1
            tind(i) = 0
            do j = 1, 12, 1
                if (tlock8(i) .eq. tglok8(j)) then
                    tind(i) = j
                    goto 30
                end if
            end do
30          continue
        end do
!
    else
        do i = 1, nblk8, 1
            tind(i) = 0
        end do
    end if
!
end subroutine
