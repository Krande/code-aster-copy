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
subroutine jm1dn2(indn, indc, nb1, nb2, xr, &
                  epais, ksi3s2, intsx, vecnph, jm1, &
                  j1dn2)
!
!
! ......................................................................
!     FONCTION :  CALCUL DU PRODUIT
!
!                 J1DN2 ( 9 ,  6 * NB1  + 3 ) =
!
!                 JTILD ( 9 , 9 ) * DNDQSI2  ( 9 ,  6 * NB1  + 3 )
!
!                 POUR LA PARTIE MATERIELLE DE LA RIGIDITE GEOMETRIQUE
!
!                 AUX POINTS D INTEGRATION REDUITE OU NORMALE
!
! ......................................................................
!
!
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/valfor.h"
    integer(kind=8) :: jn
!
    integer(kind=8) :: nb1, nb2
!
    integer(kind=8) :: indn, indc
!
    integer(kind=8) :: intsx
    integer(kind=8) :: intsx1, intsx2
!
    integer(kind=8) :: lt1, lt2
    integer(kind=8) :: i1, i2
!
    integer(kind=8) :: l1, l2, l3
    integer(kind=8) :: i3, i4, i5
!
    real(kind=8) :: vi(3)
!
    real(kind=8) :: xr(*)
!
    real(kind=8) :: epais
!
    real(kind=8) :: ksi3s2
!
    real(kind=8) :: jm1(3)
    real(kind=8) :: j1dn2(9, 51)
!
    real(kind=8) :: tmpi(3)
!
    real(kind=8) :: vecnph(9, 3)
!
!
!DEB
!
!---- INITIALISATION
!
    call r8inir(9*51, 0.d0, j1dn2, 1)
!
!---- LES ADRESSES DES FONCTIONS DE FORME ET DE LEURS DERIVEES
!     SELON INDN ( VOIR ROUTINE BTDFN )
!
    call valfor(indn, lt1, lt2, l1, l2, &
                l3)
!
!---- DECALAGE DE 8 NOEUDS DE SERENDIP
!
    intsx1 = 8*(intsx-1)
!
!---- DECALAGE DE 9 NOEUDS DE LAGRANGE
!
    intsx2 = 9*(intsx-1)
!
    i1 = lt1+intsx1
    i2 = lt2+intsx1
!
    i3 = l1+intsx2
    i4 = l2+intsx2
    i5 = l3+intsx2
!
    ASSERT((indc .eq. 1) .or. (indc .eq. 0))
!
    if (indc .eq. 1) then
!
!----- CALCUL COMPLET AVEC TERMES ROTATION
!
        do jn = 1, nb2
!
!------- PARTIE ROTATION
!
!------- REMPLISSAGE DE VI ( 3 )
!
            vi(1) = epais*ksi3s2*xr(i4+jn)
            vi(2) = epais*ksi3s2*xr(i5+jn)
            vi(3) = epais*0.5d0*xr(i3+jn)
!
!------- PRODUIT  JM1 ( 3 , 3 ) * VI ( 3 )
!
            call promat(jm1, 3, 3, 3, vi, &
                        3, 3, 1, tmpi)
!
!------- REMPLISSAGE DE J1DN2 ( 9 , 6 * NB1 + 3 )
!
!        JTILD-1 ( 9 , 9 ) * DNDQSI ( 9 , 6 * NB1 + 3 )
!
!
!
            if (jn .le. nb1) then
!
!
!
!---------- NOEUDS DE SERENDIP
!
!---------- BLOC U      0      TMPZI   -TMPYI
!
                j1dn2(1, (jn-1)*6+5) = tmpi(1)*vecnph(jn, &
                                                      3)
                j1dn2(2, (jn-1)*6+5) = tmpi(2)*vecnph(jn, &
                                                      3)
                j1dn2(3, (jn-1)*6+5) = tmpi(3)*vecnph(jn, &
                                                      3)
!
                j1dn2(1, (jn-1)*6+6) = -tmpi(1)*vecnph( &
                                       jn, 2)
                j1dn2(2, (jn-1)*6+6) = -tmpi(2)*vecnph( &
                                       jn, 2)
                j1dn2(3, (jn-1)*6+6) = -tmpi(3)*vecnph( &
                                       jn, 2)
!
!---------- BLOC V     - TMPZI 0        TMPXI
!
                j1dn2(4, (jn-1)*6+4) = -tmpi(1)*vecnph( &
                                       jn, 3)
                j1dn2(5, (jn-1)*6+4) = -tmpi(2)*vecnph( &
                                       jn, 3)
                j1dn2(6, (jn-1)*6+4) = -tmpi(3)*vecnph( &
                                       jn, 3)
!
                j1dn2(4, (jn-1)*6+6) = tmpi(1)*vecnph(jn, &
                                                      1)
                j1dn2(5, (jn-1)*6+6) = tmpi(2)*vecnph(jn, &
                                                      1)
                j1dn2(6, (jn-1)*6+6) = tmpi(3)*vecnph(jn, &
                                                      1)
!
!---------- BLOC W       TMPY  -TMPXI   0
!
                j1dn2(7, (jn-1)*6+4) = tmpi(1)*vecnph(jn, &
                                                      2)
                j1dn2(8, (jn-1)*6+4) = tmpi(2)*vecnph(jn, &
                                                      2)
                j1dn2(9, (jn-1)*6+4) = tmpi(3)*vecnph(jn, &
                                                      2)
!
                j1dn2(7, (jn-1)*6+5) = -tmpi(1)*vecnph( &
                                       jn, 1)
                j1dn2(8, (jn-1)*6+5) = -tmpi(2)*vecnph( &
                                       jn, 1)
                j1dn2(9, (jn-1)*6+5) = -tmpi(3)*vecnph( &
                                       jn, 1)
!
!
!---------- PARTIE TRANSLATION
!
!---------- REMPLISSAGE DE VI ( 3 )
!
                vi(1) = xr(i1+jn)
                vi(2) = xr(i2+jn)
                vi(3) = 0.d0
!
!---------- PRODUIT  JM1 ( 3 , 3 ) * VI ( 3 )
!
                call promat(jm1, 3, 3, 3, vi, &
                            3, 3, 1, tmpi)
!
!---------- BLOC U
!
                j1dn2(1, (jn-1)*6+1) = tmpi(1)
                j1dn2(2, (jn-1)*6+1) = tmpi(2)
                j1dn2(3, (jn-1)*6+1) = tmpi(3)
!
!---------- BLOC V
!
                j1dn2(4, (jn-1)*6+2) = tmpi(1)
                j1dn2(5, (jn-1)*6+2) = tmpi(2)
                j1dn2(6, (jn-1)*6+2) = tmpi(3)
!
!---------- BLOC W
!
                j1dn2(7, (jn-1)*6+3) = tmpi(1)
                j1dn2(8, (jn-1)*6+3) = tmpi(2)
                j1dn2(9, (jn-1)*6+3) = tmpi(3)
!
            else
!
!------- SUPERNOEUD
!
!
!
!
!---------- BLOC U      0      TMPZI   -TMPYI
!
                j1dn2(1, nb1*6+2) = tmpi(1)*vecnph(jn, 3 &
                                                   )
                j1dn2(2, nb1*6+2) = tmpi(2)*vecnph(jn, 3 &
                                                   )
                j1dn2(3, nb1*6+2) = tmpi(3)*vecnph(jn, 3 &
                                                   )
!
                j1dn2(1, nb1*6+3) = -tmpi(1)*vecnph(jn, &
                                                    2)
                j1dn2(2, nb1*6+3) = -tmpi(2)*vecnph(jn, &
                                                    2)
                j1dn2(3, nb1*6+3) = -tmpi(3)*vecnph(jn, &
                                                    2)
!
!---------- BLOC V     - TMPZI 0        TMPXI
!
                j1dn2(4, nb1*6+1) = -tmpi(1)*vecnph(jn, &
                                                    3)
                j1dn2(5, nb1*6+1) = -tmpi(2)*vecnph(jn, &
                                                    3)
                j1dn2(6, nb1*6+1) = -tmpi(3)*vecnph(jn, &
                                                    3)
!
                j1dn2(4, nb1*6+3) = tmpi(1)*vecnph(jn, 1 &
                                                   )
                j1dn2(5, nb1*6+3) = tmpi(2)*vecnph(jn, 1 &
                                                   )
                j1dn2(6, nb1*6+3) = tmpi(3)*vecnph(jn, 1 &
                                                   )
!
!---------- BLOC W       TMPY  -TMPXI   0
!
                j1dn2(7, nb1*6+1) = tmpi(1)*vecnph(jn, 2 &
                                                   )
                j1dn2(8, nb1*6+1) = tmpi(2)*vecnph(jn, 2 &
                                                   )
                j1dn2(9, nb1*6+1) = tmpi(3)*vecnph(jn, 2 &
                                                   )
!
                j1dn2(7, nb1*6+2) = -tmpi(1)*vecnph(jn, &
                                                    1)
                j1dn2(8, nb1*6+2) = -tmpi(2)*vecnph(jn, &
                                                    1)
                j1dn2(9, nb1*6+2) = -tmpi(3)*vecnph(jn, &
                                                    1)
!
            end if
!
        end do
!
!
!
!
    else if (indc .eq. 0) then
!
!
!
!
!----- CALCUL INCOMPLET SANS TERMES ROTATION
!
!
        do jn = 1, nb1
!
!---------- NOEUDS DE SERENDIP
!
!---------- PARTIE TRANSLATION
!
!---------- REMPLISSAGE DE VI ( 3 )
!
            vi(1) = xr(i1+jn)
            vi(2) = xr(i2+jn)
            vi(3) = 0.d0
!
!---------- PRODUIT  JM1 ( 3 , 3 ) * VI ( 3 )
!
            call promat(jm1, 3, 3, 3, vi, &
                        3, 3, 1, tmpi)
!
!---------- BLOC U
!
            j1dn2(1, (jn-1)*6+1) = tmpi(1)
            j1dn2(2, (jn-1)*6+1) = tmpi(2)
            j1dn2(3, (jn-1)*6+1) = tmpi(3)
!
!---------- BLOC V
!
            j1dn2(4, (jn-1)*6+2) = tmpi(1)
            j1dn2(5, (jn-1)*6+2) = tmpi(2)
            j1dn2(6, (jn-1)*6+2) = tmpi(3)
!
!---------- BLOC W
!
            j1dn2(7, (jn-1)*6+3) = tmpi(1)
            j1dn2(8, (jn-1)*6+3) = tmpi(2)
            j1dn2(9, (jn-1)*6+3) = tmpi(3)
!
        end do
!
!
!
    end if
!
!
!
!
!FIN
!
end subroutine
