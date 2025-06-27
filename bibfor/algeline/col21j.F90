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
subroutine col21j(fronti, frontj, frn, j, l, &
                  n, n1, t1, t2)
! person_in_charge: olivier.boiteau at edf.fr
!     VERSION MODIFIEE POUR L' APPEL A DGEMV (PRODUITS MATRICE-VECTEUR)
!     LE STOCKAGE DES COLONNES DE LA FACTORISEE EST MODIFIE
    implicit none
    integer(kind=8) :: n, n1
    real(kind=8) :: fronti(*), frontj(*), t1(n), t2(n), frn(*)
!
    integer(kind=8) :: i, j, l, ll, ic1, ic2, id1, id2, jd1, jd
!-----------------------------------------------------------------------
    integer(kind=8) :: k
!-----------------------------------------------------------------------
    ic1 = 3
    ic2 = ic1+n
    jd1 = 0
    ll = l
    do k = 1, n1
        id1 = ic1
        id2 = ic2
        do i = 1, ll
            jd1 = jd1+1
            frontj(jd1) = frontj(jd1)-t1(k)*fronti(id1)-t2(k)*fronti(id2)
            id1 = id1+1
            id2 = id2+1
        end do
        ll = ll-1
        ic1 = ic1+1
        ic2 = ic2+1
        jd1 = jd1+2*j+k
    end do
    jd1 = 0
    do k = n1+1, l
        id1 = ic1
        id2 = ic2
        jd = jd1
        do i = 1, ll
            jd = jd+1
            frn(jd) = frn(jd)-t1(k)*fronti(id1)
            id1 = id1+1
        end do
!
        do i = 1, ll
            jd1 = jd1+1
            frn(jd1) = frn(jd1)-t2(k)*fronti(id2)
            id2 = id2+1
        end do
        ll = ll-1
        ic1 = ic1+1
        ic2 = ic2+1
    end do
end subroutine
