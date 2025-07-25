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
subroutine btdfn(ind, nb1, nb2, ksi3s2, intsn, &
                 xr, epais, vectpt, hsj1fx, btdf)
!
    implicit none
    integer(kind=8) :: nb1, nb2, intsn
    real(kind=8) :: xr(*), epais, vectpt(9, 2, 3)
    real(kind=8) :: hsj1fx(3, 9), btdf(3, 42)
    real(kind=8) :: dnsdsf(9, 42)
    real(kind=8) :: ksi3s2
    common/dnsf/dnsdsf
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i3, i4, i5, ind, intsn1, j
    integer(kind=8) :: j1, jb, k, l1, l2, l3
!-----------------------------------------------------------------------
    if (ind .eq. 1) then
        l1 = 459
        l2 = 540
        l3 = 621
    else if (ind .eq. 0) then
        l1 = 351
        l2 = 387
        l3 = 423
    end if
!
    do j = 1, 5*nb1+2
        do i = 1, 9
            dnsdsf(i, j) = 0.d0
        end do
        do i = 1, 3
            btdf(i, j) = 0.d0
        end do
    end do
!
!
    intsn1 = 9*(intsn-1)
!
!                         DN
!     CONSTRUCTION DE   ------    AUX PTS DE GAUSS NORMAL
!                        DQSI  F
!
    i3 = l1+intsn1
    i4 = l2+intsn1
    i5 = l3+intsn1
    do j = 1, nb1
        j1 = 5*(j-1)
        dnsdsf(1, j1+4) = -ksi3s2*xr(i4+j)*epais*vectpt(j, 2, 1)
        dnsdsf(1, j1+5) = ksi3s2*xr(i4+j)*epais*vectpt(j, 1, 1)
!
        dnsdsf(2, j1+4) = -ksi3s2*xr(i5+j)*epais*vectpt(j, 2, 1)
        dnsdsf(2, j1+5) = ksi3s2*xr(i5+j)*epais*vectpt(j, 1, 1)
!
        dnsdsf(3, j1+4) = -xr(i3+j)/2*epais*vectpt(j, 2, 1)
        dnsdsf(3, j1+5) = xr(i3+j)/2*epais*vectpt(j, 1, 1)
!
        dnsdsf(4, j1+4) = -ksi3s2*xr(i4+j)*epais*vectpt(j, 2, 2)
        dnsdsf(4, j1+5) = ksi3s2*xr(i4+j)*epais*vectpt(j, 1, 2)
!
        dnsdsf(5, j1+4) = -ksi3s2*xr(i5+j)*epais*vectpt(j, 2, 2)
        dnsdsf(5, j1+5) = ksi3s2*xr(i5+j)*epais*vectpt(j, 1, 2)
!
        dnsdsf(6, j1+4) = -xr(i3+j)/2*epais*vectpt(j, 2, 2)
        dnsdsf(6, j1+5) = xr(i3+j)/2*epais*vectpt(j, 1, 2)
!
        dnsdsf(7, j1+4) = -ksi3s2*xr(i4+j)*epais*vectpt(j, 2, 3)
        dnsdsf(7, j1+5) = ksi3s2*xr(i4+j)*epais*vectpt(j, 1, 3)
!
        dnsdsf(8, j1+4) = -ksi3s2*xr(i5+j)*epais*vectpt(j, 2, 3)
        dnsdsf(8, j1+5) = ksi3s2*xr(i5+j)*epais*vectpt(j, 1, 3)
!
        dnsdsf(9, j1+4) = -xr(i3+j)/2*epais*vectpt(j, 2, 3)
        dnsdsf(9, j1+5) = xr(i3+j)/2*epais*vectpt(j, 1, 3)
    end do
!
    dnsdsf(1, 5*nb1+1) = -ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 2, 1)
    dnsdsf(1, 5*nb1+2) = ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 1, 1)
!
    dnsdsf(2, 5*nb1+1) = -ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 2, 1)
    dnsdsf(2, 5*nb1+2) = ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 1, 1)
!
    dnsdsf(3, 5*nb1+1) = -xr(i3+nb2)/2*epais*vectpt(nb2, 2, 1)
    dnsdsf(3, 5*nb1+2) = xr(i3+nb2)/2*epais*vectpt(nb2, 1, 1)
!
    dnsdsf(4, 5*nb1+1) = -ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 2, 2)
    dnsdsf(4, 5*nb1+2) = ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 1, 2)
!
    dnsdsf(5, 5*nb1+1) = -ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 2, 2)
    dnsdsf(5, 5*nb1+2) = ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 1, 2)
!
    dnsdsf(6, 5*nb1+1) = -xr(i3+nb2)/2*epais*vectpt(nb2, 2, 2)
    dnsdsf(6, 5*nb1+2) = xr(i3+nb2)/2*epais*vectpt(nb2, 1, 2)
!
    dnsdsf(7, 5*nb1+1) = -ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 2, 3)
    dnsdsf(7, 5*nb1+2) = ksi3s2*xr(i4+nb2)*epais*vectpt(nb2, 1, 3)
!
    dnsdsf(8, 5*nb1+1) = -ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 2, 3)
    dnsdsf(8, 5*nb1+2) = ksi3s2*xr(i5+nb2)*epais*vectpt(nb2, 1, 3)
!
    dnsdsf(9, 5*nb1+1) = -xr(i3+nb2)/2*epais*vectpt(nb2, 2, 3)
    dnsdsf(9, 5*nb1+2) = xr(i3+nb2)/2*epais*vectpt(nb2, 1, 3)
!
!     CONSTRUCTION DE BTILDF = HFM * S * JTILD-1 * DNSDSF  : (3,5*NB1+2)
!
    do i = 1, 3
        do jb = 1, nb1
            do j = 4, 5
                j1 = j+5*(jb-1)
                btdf(i, j1) = 0.d0
                do k = 1, 9
                    btdf(i, j1) = btdf(i, j1)+hsj1fx(i, k)*dnsdsf(k, j1)
                end do
            end do
        end do
!
        do j = 1, 2
            j1 = j+5*nb1
            btdf(i, j1) = 0.d0
            do k = 1, 9
                btdf(i, j1) = btdf(i, j1)+hsj1fx(i, k)*dnsdsf(k, j1)
            end do
        end do
    end do
!
end subroutine
