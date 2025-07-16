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
subroutine vpzbal(mat, neq, mxeq, d, k, &
                  l)
    implicit none
#include "asterc/r8baem.h"
    integer(kind=8) :: neq, mxeq, k, l
    real(kind=8) :: mat(mxeq, 1), d(1)
!     REDUCTION DE LA NORME DE LA MATRICE PAR LA TRANSFORMATION DE
!     SIMILITUDE STOCKEE DANS "D"
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE 320
!     ------------------------------------------------------------------
    integer(kind=8) :: l1, k1, j, i, ll, noconv
    real(kind=8) :: b, b2, r, c, f, g, s
!     ------------------------------------------------------------------
!     --- RECUPERATION DE LA BASE DE NUMEROTATION DE LA MACHINE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    b = r8baem()
    b2 = b*b
!
!     ---RECHERCHE DES VALEURS PROPRES ISOLEES (ON LES METS A LA FIN)---
!     --- CAS DES LIGNES ---
    l1 = 1
    k1 = neq
5   continue
    do j = k1, 1, -1
        r = -abs(mat(j, j))
        do i = 1, k1
            r = r+abs(mat(j, i))
        end do
        if (r .eq. 0.d0) then
            d(k1) = j
            if (j .ne. k1) then
                do i = 1, k1
                    f = mat(i, j)
                    mat(i, j) = mat(i, k1)
                    mat(i, k1) = f
                end do
                do i = l1, neq
                    f = mat(j, i)
                    mat(j, i) = mat(k1, i)
                    mat(k1, i) = f
                end do
            end if
            k1 = k1-1
            goto 5
        end if
    end do
!
!     ---RECHERCHE DES VALEURS PROPRES ISOLEES (ON LES METS A GAUCHE)---
!     --- CAS DES COLONNES -
35  continue
    ll = l1
    do j = ll, k1
        c = -abs(mat(j, j))
        do i = l1, k1
            c = c+abs(mat(i, j))
        end do
        if (c .eq. 0.d0) then
            d(l1) = j
            if (j .ne. l1) then
                do i = 1, k1
                    f = mat(i, j)
                    mat(i, j) = mat(i, l1)
                    mat(i, l1) = f
                end do
                do i = l1, neq
                    f = mat(j, i)
                    mat(j, i) = mat(l1, i)
                    mat(l1, i) = f
                end do
            end if
            l1 = l1+1
            goto 35
        end if
    end do
!
!     EQUILIBRER LA SOUS-MATRICE DE LA LIGNES L1 A K1
    k = l1
    l = k1
    do i = l1, k1
        d(i) = 1.d0
    end do
75  continue
    noconv = 0
    do i = l1, k1
        c = -abs(mat(i, i))
        r = c
        do j = l1, k1
            c = c+abs(mat(j, i))
            r = r+abs(mat(i, j))
        end do
        g = r/b
        f = 1.d0
        s = c+r
85      continue
        if (c .lt. g) then
            f = f*b
            c = c*b2
            goto 85
        end if
        g = r*b
95      continue
        if (c .ge. g) then
            f = f/b
            c = c/b2
            goto 95
        end if
!
!        --- EQUILIBRAGE ---
        if ((c+r)/f .lt. 0.95d0*s) then
            g = 1.d0/f
            d(i) = d(i)*f
            noconv = 1
            do j = l1, neq
                mat(i, j) = mat(i, j)*g
            end do
            do j = 1, k1
                mat(j, i) = mat(j, i)*f
            end do
        end if
    end do
    if (noconv .eq. 1) goto 75
end subroutine
