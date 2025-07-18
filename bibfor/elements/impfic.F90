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
subroutine impfic(vale, nomnoe, rcmp, unit, lxfem)
    implicit none
#include "asterf_types.h"
! ......................................................................
!
!     - FONCTION REALISEE:  IMPRESSION DES FACTEURS D'INTENSITE
!                           DE CONTRAINTES K1 ET K2, DU TAUX DE
!                           RESTITUTION D'ENERGIE G (VALEUR CLASSIQUE G
!                           ET PAR LA FORMULE D'IRWIN )
!
!  ENTREE
!
!    VALE       --> ADRESSE DES VALEURS DES FIC
!    NOMNOE     --> ADRESSE DU NOM DU NOEUD DE FOND DE FISSURE
!    RCMP       --> COORDONNEES DU NOEUD DE FOND DE FISSURE
!                     ET DE LA NORMALE A LA FISSURE
!    LXFEM       --> FALSE SI ON EST DANS UN CALCUL CLASSIQUE (FEM)
!                    TRUE SI ON EST DANS UN CALCUL X-FEM
! ......................................................................
!
    integer(kind=8) :: i, k1phi, k2phi, gphi, unit
    real(kind=8) :: g, fic1, fic2, k1, k2, girwin, rcmp(6), vale(5)
    real(kind=8) :: k11(10), k21(10), k12(10), k22(10), k1max, k1min
    real(kind=8) :: k2min, k2sup, gmax, k1dev, k2dev, gdev, fic1d, fic2d
    character(len=8) :: nomnoe
    aster_logical :: lxfem
! ......................................................................
!
    g = vale(1)
    fic1 = vale(2)
    fic2 = vale(3)
    k1 = vale(4)
    k2 = vale(5)
    girwin = fic1*fic1+fic2*fic2
!
! CALCUL DE L'ANGLE DE PROPAGATION DE LA FISSURE SELON 3 CRITERES:
!  K1MAX, K2MIN, GMAX (FORMULES ISSUES DE AMNESTOY-BUI)
!
    k11(1) = 1.d0
    k11(2) = 0.9886d0
    k11(3) = 0.9552d0
    k11(4) = 0.9018d0
    k11(5) = 0.8314d0
    k11(6) = 0.7479d0
    k11(7) = 0.6559d0
    k11(8) = 0.5598d0
    k11(9) = 0.4640d0
    k11(10) = 0.3722d0
!
    k21(1) = 0
    k21(2) = 0.0864d0
    k21(3) = 0.1680d0
    k21(4) = 0.2403d0
    k21(5) = 0.2995d0
    k21(6) = 0.3431d0
    k21(7) = 0.3696d0
    k21(8) = 0.3788d0
    k21(9) = 0.3718d0
    k21(10) = 0.3507d0
!
    k12(1) = 0.d0
    k12(2) = -0.2597d0
    k12(3) = -0.5068d0
    k12(4) = -0.7298d0
    k12(5) = -0.9189d0
    k12(6) = -1.0665d0
    k12(7) = -1.1681d0
    k12(8) = -1.2220d0
    k12(9) = -1.2293d0
    k12(10) = -1.1936d0
!
    k22(1) = 1.d0
    k22(2) = 0.9764d0
    k22(3) = 0.9071d0
    k22(4) = 0.7972d0
    k22(5) = 0.6540d0
    k22(6) = 0.4872d0
    k22(7) = 0.3077d0
    k22(8) = 0.1266d0
    k22(9) = -0.0453d0
    k22(10) = -0.1988d0
!
    k1max = 0.d0
    k1min = 0.d0
    k2min = -abs(k2)
    k2sup = +abs(k2)
    gmax = 0.d0
    k1phi = 0
    k2phi = 0
    gphi = 0
!
    do i = 1, 10
        k1dev = k11(i)*k1+k12(i)*k2
        k2dev = k21(i)*k1+k22(i)*k2
        fic1d = k11(i)*fic1+k12(i)*fic2
        fic2d = k21(i)*fic1+k22(i)*fic2
        gdev = fic1d*fic1d+fic2d*fic2d
        if ((k1dev .gt. k1max) .or. (k1dev .lt. k1min)) then
            k1phi = 10*(i-1)
            k1min = -abs(k1dev)
            k1max = abs(k1dev)
        end if
        if ((k2dev .le. k2sup) .and. (k2dev .ge. k2min)) then
            k2phi = 10*(i-1)
            k2min = -abs(k2dev)
            k2sup = abs(k2dev)
        end if
        if (gdev .gt. gmax) then
            gphi = 10*(i-1)
            gmax = gdev
        end if
    end do
    do i = 2, 10
        k1dev = k11(i)*k1-k12(i)*k2
        k2dev = -k21(i)*k1+k22(i)*k2
        fic1d = k11(i)*fic1-k12(i)*fic2
        fic2d = -k21(i)*fic1+k22(i)*fic2
        gdev = fic1d*fic1d+fic2d*fic2d
        if ((k1dev .gt. k1max) .or. (k1dev .lt. k1min)) then
            k1phi = -10*(i-1)
            k1min = -abs(k1dev)
            k1max = abs(k1dev)
        end if
        if ((k2dev .le. k2sup) .and. (k2dev .ge. k2min)) then
            k2phi = -10*(i-1)
            k2min = -abs(k2dev)
            k2sup = abs(k2dev)
        end if
        if (gdev .gt. gmax) then
            gphi = -10*(i-1)
            gmax = gdev
        end if
    end do
!
    write (unit, *)
    write (unit, 555)
    write (unit, *)
!
    if (.not. lxfem) then
        write (unit, *) 'NOEUD DE FOND DE FISSURE : ', nomnoe
    end if
    write (unit, *)
!
    write (unit, *) 'COORDONNEES DU NOEUD DE FOND DE FISSURE : ',&
     &              rcmp(1), ' ', rcmp(2)
    write (unit, *)
    write (unit, *) 'COORDONNEES DE LA NORMALE A LA FISSURE :  ',&
     &              rcmp(5), ' ', rcmp(6)
    write (unit, *)
!
    write (unit, *) '       K1                K2         '//&
     &               '     G (IRWIN)'
    write (unit, *)
    write (unit, 999) k1, k2, girwin
    write (unit, *)
    write (unit, 777) 'TAUX DE RESTITUTION D''ENERGIE G : ', g
    write (unit, *)
    write (unit, *) 'DIRECTION DE LA DEVIATION DE LA FISSURE '&
     &             //'(EN DEGRES): '
    write (unit, *)
    write (unit, 666) 'SELON LE CRITERE K1 MAXIMUM : ', k1phi,&
     &              ' AVEC K1MAX = ', k1max
    write (unit, 666) 'SELON LE CRITERE K2 NUL     : ', k2phi,&
     &              ' AVEC K2NUL = ', k2sup
    write (unit, 666) 'SELON LE CRITERE G MAXIMUM  : ', gphi,&
     &              ' AVEC GMAX  = ', gmax
    write (unit, 555)
    write (unit, *)
!
555 format(60('*'))
666 format(a, i3, a, 1pd12.5)
777 format(a, 1pd12.5)
999 format(3(1pd12.5, 8x))
!
end subroutine
