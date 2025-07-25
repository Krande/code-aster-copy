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

subroutine dkqshp(qsi, eta, caraq4, jacob, shpN, shpr1, shpr2)
    implicit none

#include "asterf_types.h"
    real(kind=8) :: qsi, eta, caraq4(*), jacob(*)
    real(kind=8) :: shpN(3, 4), shpr1(3, 4), shpr2(3, 4)

!
!     ------------------------------------------------------------------
!     MATRICE SHP(3,4) DES FONCTIONS DE BASE
!     AU POINT QSI ETA POUR ELEMENTS DKQ ET DSQ
!     ------------------------------------------------------------------

    integer(kind=8) :: i, j, l
    real(kind=8) :: vj11, vj12, vj21, vj22
    real(kind=8) :: s2, t2, shp(3, 4), temp, n1(4), n2(4)
!     ------------------------------------------------------------------
!   Jacobian
!
    vj11 = jacob(1)
    vj12 = jacob(2)
    vj21 = jacob(3)
    vj22 = jacob(4)
!
!   Form bilinear shape functions Ni
!

!
!   Ni
!
    shpN(3, 1) = 0.25d0*(1.d0-eta)*(1.d0-qsi)
    shpN(3, 2) = 0.25d0*(1.d0-eta)*(1.d0+qsi)
    shpN(3, 3) = 0.25d0*(1.d0+eta)*(1.d0+qsi)
    shpN(3, 4) = 0.25d0*(1.d0+eta)*(1.d0-qsi)
!
!   dNi/d_qsi
!
    shpN(1, 1) = -0.25d0*(1.d0-eta)
    shpN(1, 2) = 0.25d0*(1.d0-eta)
    shpN(1, 3) = 0.25d0*(1.d0+eta)
    shpN(1, 4) = -0.25d0*(1.d0+eta)
!
!   dNi/d_eta
!
    shpN(2, 1) = -0.25d0*(1.d0-qsi)
    shpN(2, 2) = -0.25d0*(1.d0+qsi)
    shpN(2, 3) = 0.25d0*(1.d0+qsi)
    shpN(2, 4) = 0.25d0*(1.d0-qsi)

!
!   Form global derivatives
!
    do i = 1, 4
        temp = shpN(1, i)*vj11+shpN(2, i)*vj12
        shpN(2, i) = shpN(1, i)*vj21+shpN(2, i)*vj22
        shpN(1, i) = temp
    end do

!
!   Form quadratic shape functions
!
    s2 = (1.d0-qsi*qsi)/2.0d0
    t2 = (1.d0-eta*eta)/2.0d0
!
!   Si
!
    shp(3, 1) = s2*(1.d0-eta)
    shp(3, 2) = t2*(1.d0+qsi)
    shp(3, 3) = s2*(1.d0+eta)
    shp(3, 4) = t2*(1.d0-qsi)
!
!   dSi/d_qsi
!
    shp(1, 1) = -qsi*(1.d0-eta)
    shp(1, 2) = t2
    shp(1, 3) = -qsi*(1.d0+eta)
    shp(1, 4) = -t2
!
!   dSi/d_eta
!
    shp(2, 1) = -s2
    shp(2, 2) = -eta*(1.d0+qsi)
    shp(2, 3) = s2
    shp(2, 4) = -eta*(1.d0-qsi)
!
!   Form global derivatives
!
    do i = 1, 4
        temp = shp(1, i)*vj11+shp(2, i)*vj12
        shp(2, i) = shp(1, i)*vj21+shp(2, i)*vj22
        shp(1, i) = temp
    end do

!
!   Form rotational shape functions
!
!

    n2(1) = -caraq4(13)*caraq4(9)*0.125d0
    n2(2) = -caraq4(14)*caraq4(10)*0.125d0
    n2(3) = -caraq4(15)*caraq4(11)*0.125d0
    n2(4) = -caraq4(16)*caraq4(12)*0.125d0

    n1(1) = +caraq4(17)*caraq4(9)*0.125d0
    n1(2) = +caraq4(18)*caraq4(10)*0.125d0
    n1(3) = +caraq4(19)*caraq4(11)*0.125d0
    n1(4) = +caraq4(20)*caraq4(12)*0.125d0

    do l = 1, 3
        j = 4
        do i = 1, 4

            shpr1(l, i) = shp(l, j)*n1(j)-shp(l, i)*n1(i)
            shpr2(l, i) = shp(l, j)*n2(j)-shp(l, i)*n2(i)
            j = i
        end do
    end do
end subroutine
