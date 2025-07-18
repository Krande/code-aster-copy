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
subroutine angco4(coor, zk1, izk, icoude, zk2, &
                  rayon, theta, angl1, angl2, angl3, &
                  angl4, pgl1, pgl2, pgl3, pgl4, &
                  omega, dn1n2, epsi, crit)
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/angcou.h"
#include "asterfort/assert.h"
    real(kind=8) :: coor(*), rayon, theta, epsi
    real(kind=8) :: coor3(12)
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), valtes
!     ON POURRAIT SE PASSER DE PGLI
    real(kind=8) :: angl1(3), angl2(3), angl3(3), rayon1, omega1
    real(kind=8) :: rayon2, theta1, theta2, angl4(3), pgl4(3, 3), test
    real(kind=8) :: zk1(3), zk2(3), zkini(3), zk4(3), zk3(3), epsi2
    real(kind=8) :: omega, dn1n2, omega2, coo1(3), coo2(3), dn1n4, dn3n2
    character(len=8) :: crit
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DE LA GEOMETRIE COUDE
!                          TUYAU
!    - ARGUMENTS
!        DONNEES:      COOR       -->  CORDONNEES DE SNOEUDS
!         SORTIE:      ICOUDE      =0 DROIT =1 COUDE
!         SORTIE:      ANGL1,2,3
! ......................................................................
!
    integer(kind=8) :: icoude, i, izk, icoud1, icoud2
!
!     POUR VERIFICATIONS (PAS TRES EXIGEANTES) SUR LA GEOMETRIE
!     SINON, IL FAUDRAIT INTRODUIRE UN AUTRE MOT CLE PRECISON2
!     DIFFERENT DE PRECISION
!
    epsi2 = 1.d-4
!
!     POINTS 1 4 3
!
    do i = 1, 3
        coor3(i) = coor(i)
        coor3(3+i) = coor(9+i)
        coor3(6+i) = coor(6+i)
        zkini(i) = zk1(i)
    end do
!
!
    call angcou(coor3, zkini, izk, icoud1, zk4, &
                rayon1, theta1, angl1, angl4, angl3, &
                pgl1, pgl4, pgl3, omega1, dn1n4, &
                epsi, crit, zk3)
!
!     POINTS 3 2 4
!
    do i = 1, 3
        coor3(i) = coor(6+i)
        coor3(3+i) = coor(3+i)
        coor3(6+i) = coor(9+i)
        zkini(i) = zk3(i)
    end do
!
    call angcou(coor3, zkini, izk, icoud2, zk2, &
                rayon2, theta2, angl3, angl2, angl4, &
                pgl3, pgl2, pgl4, omega2, dn3n2, &
                epsi, crit, zk4)
!
    do i = 1, 3
        coo1(i) = coor(i)
        coo2(i) = coor(3+i)
    end do
    dn1n2 = sqrt((coo1(1)-coo2(1))**2+(coo1(2)-coo2(2))**2+(coo1(3)-coo2(3))**2)
!
    if (icoud2 .ne. icoud1) then
        ASSERT(.false.)
    else
        icoude = icoud2
        rayon = rayon2
    end if
!
    valtes = abs(rayon1)
    if (crit .eq. 'RELATIF') then
        if (valtes .lt. r8prem()) then
            test = r8prem()
        else
            test = epsi2*valtes
        end if
    else if (crit .eq. 'ABSOLU') then
        test = epsi2
    end if
!
    if (abs(rayon2-rayon1) .gt. test) then
        ASSERT(.false.)
    else
        rayon = rayon2
    end if
!
    valtes = abs(theta1)
    if (crit .eq. 'RELATIF') then
        if (valtes .lt. r8prem()) then
            test = r8prem()
        else
            test = epsi2*valtes
        end if
    else if (crit .eq. 'ABSOLU') then
        test = epsi2
    end if
!
    if (abs(theta2-theta1) .gt. test) then
        ASSERT(.false.)
    else
        theta = 1.5d0*theta1
    end if
!
    valtes = abs(omega1)
    if (crit .eq. 'RELATIF') then
        if (valtes .lt. r8prem()) then
            test = r8prem()
        else
            test = epsi2*valtes
        end if
    else if (crit .eq. 'ABSOLU') then
        test = epsi2
    end if
!
    if (abs(omega2-omega1) .gt. test) then
        ASSERT(.false.)
    else
        omega = omega1
    end if
end subroutine
