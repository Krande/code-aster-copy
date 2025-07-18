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
subroutine vectgt(ind, nb1, xi, ksi3s2, intsx, &
                  zr, epais, vectn, vectg, vectt)
!
    implicit none
    integer(kind=8) :: ind, nb1, intsx
    real(kind=8) :: xi(3, *), zr(*), epais, vectn(9, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: ksi3s2
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, intsx1, j, k, l1, l2
    integer(kind=8) :: l3
    real(kind=8) :: rnorm
!-----------------------------------------------------------------------
    if (ind .eq. 0) then
!
!     CALCULS AUX PTS D'INTEGRATION REDUITE
!
        l1 = 12
        l2 = 44
        l3 = 76
    else if (ind .eq. 1) then
!
!     CALCULS AUX PTS D'INTEGRATION NORMALE
!
        l1 = 135
        l2 = 207
        l3 = 279
!
    end if
!
!     CONSTRUCTION DU VECTEUR N AUX X PTS DE GAUSS. X = REDUIT OU NORMAL
!     (STOCKE DANS VECTT)
!
    intsx1 = 8*(intsx-1)
    i1 = l1+intsx1
    do k = 1, 3
        vectt(3, k) = 0
        do j = 1, nb1
            vectt(3, k) = vectt(3, k)+zr(i1+j)*vectn(j, k)
        end do
    end do
!
!     CONSTRUCTION DES VECTEURS GA AUX X PTS DE GAUSS
!
    i1 = l2+intsx1
    i2 = l3+intsx1
    do k = 1, 3
        vectg(1, k) = 0.d0
        vectg(2, k) = 0.d0
        do j = 1, nb1
            vectg(1, k) = vectg(1, k)+zr(i1+j)*(xi(k, j)+ksi3s2*epais* &
                                                vectn(j, k))
            vectg(2, k) = vectg(2, k)+zr(i2+j)*(xi(k, j)+ksi3s2*epais* &
                                                vectn(j, k))
        end do
    end do
!
!     CONSTRUCTION DES VECTEURS TA AUX X PTS DE GAUSS (T3=N)
!
    rnorm = sqrt(vectg(1, 1)*vectg(1, 1)&
     &             +vectg(1, 2)*vectg(1, 2)&
     &             +vectg(1, 3)*vectg(1, 3))
!
    do k = 1, 3
        vectt(1, k) = vectg(1, k)/rnorm
    end do
!
    vectt(2, 1) = vectt(3, 2)*vectt(1, 3)-vectt(3, 3)*vectt(1, 2)
    vectt(2, 2) = vectt(3, 3)*vectt(1, 1)-vectt(3, 1)*vectt(1, 3)
    vectt(2, 3) = vectt(3, 1)*vectt(1, 2)-vectt(3, 2)*vectt(1, 1)
!
end subroutine
