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
subroutine vpzhes(mat, k, l, neq, mxeq, &
                  d)
    implicit none
    integer(kind=8) :: k, l, neq, mxeq
    real(kind=8) :: mat(mxeq, neq), d(neq)
!     MISE SOUS FORME DE HESSENBERG (FORME SUPERIEURE)
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE 342
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: m, i, j
    real(kind=8) :: f, g, h, scale, zero
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    zero = 0.d0
    do m = k+1, l-1
        h = zero
        d(m) = zero
        scale = zero
!
!        --- MISE A L'ECHELLE DE LA COLONNE ---
        do i = m, l
            scale = scale+abs(mat(i, m-1))
        end do
        if (scale .eq. zero) goto 45
        do i = l, m, -1
            d(i) = mat(i, m-1)/scale
            h = h+d(i)*d(i)
        end do
        g = -sign(sqrt(h), d(m))
        h = h-d(m)*g
        d(m) = d(m)-g
!
!        --- FORMATION DE  (I-(U*UT)/H) * MAT  ---
        do j = m, neq
            f = zero
            do i = l, m, -1
                f = f+d(i)*mat(i, j)
            end do
            f = f/h
            do i = m, l
                mat(i, j) = mat(i, j)-f*d(i)
            end do
        end do
!
!        --- FORMATION DE (I-(U*UT)/H)*MAT*(I-(U*UT)/H)  ---
        do i = 1, l
            f = zero
            do j = l, m, -1
                f = f+d(j)*mat(i, j)
            end do
            f = f/h
            do j = m, l
                mat(i, j) = mat(i, j)-f*d(j)
            end do
        end do
        d(m) = scale*d(m)
        mat(m, m-1) = scale*g
45      continue
    end do
end subroutine
