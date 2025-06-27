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
subroutine vpzrbk(z, h, d, mm, izh, &
                  k, l)
    implicit none
    integer(kind=8) :: mm, izh, k, l
    real(kind=8) :: z(izh, 1), h(izh, 1), d(1)
!     TRANSFORMATION ARRIERE POUR OBTENIR LES VECTEURS (ROUTINE ORTBAK)
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE 350
!     ------------------------------------------------------------------
    integer(kind=8) :: m, ma, i, j
    real(kind=8) :: g, zero
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    zero = 0.d0
    do m = l-2, k, -1
        ma = m+1
        if (h(ma, m) .ne. zero) then
            do i = m+2, l
                d(i) = h(i, m)
            end do
            if (ma .le. l) then
                do j = 1, mm
                    g = zero
                    do i = ma, l
                        g = g+d(i)*z(i, j)
                    end do
!
                    g = (g/d(ma))/h(ma, m)
                    do i = ma, l
                        z(i, j) = z(i, j)+g*d(i)
                    end do
                end do
            end if
        end if
    end do
end subroutine
