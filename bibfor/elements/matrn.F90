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
subroutine matrn(nb1, nb2, xr, ksi3s2, epais, &
                 intsn, vectn, matn)
!
!
! ......................................................................
!     FONCTION :  CALCUL DE
!
!                 MATN ( 3 ,  6 * NB1  + 3 ) =
!
!                 AUX POINTS D INTEGRATION NORMALE
!
! ......................................................................
!
!
!
    implicit none
!
#include "asterfort/antisy.h"
#include "asterfort/r8inir.h"
    real(kind=8) :: matn(3, 51)
    real(kind=8) :: xr(*)
    real(kind=8) :: vectn(9, 3)
!
    integer(kind=8) :: jn
!
    integer(kind=8) :: intsn
!
    integer(kind=8) :: nb1, nb2
    integer(kind=8) :: ii, jj
!
    real(kind=8) :: vecnj(3), antnj(3, 3)
    real(kind=8) :: ksi3s2
    real(kind=8) :: epais
!
!DEB
!
!---- INITIALISATION
!
    call r8inir(3*51, 0.d0, matn, 1)
!
!---- LES ADRESSES DES FONCTIONS DE FORME ET DE LEURS DERIVEES
!     DECALAGE DE 8 NOEUDS DE SERENDIP ET 9 NOEUDS DE LAGRANGE
!
!
    do jn = 1, nb2
!
!------- NORMALE ET ANTISYM AU NOEUD JN
!
        do ii = 1, 3
            vecnj(ii) = vectn(jn, ii)
        end do
!
        call antisy(vecnj, 1.d0, antnj)
!
        if (jn .le. nb1) then
!
!---------- NOEUDS DE SERENDIP
!
!---------- PARTIE TRANSLATION
!
            do jj = 1, 3
                matn(jj, (jn-1)*6+jj) = xr(135+8*( &
                                           intsn-1)+jn)
            end do
!
!---------- PARTIE ROTATION
!
            do jj = 1, 3
                do ii = 1, 3
                    matn(ii, (jn-1)*6+jj+3) = -epais* &
                                              ksi3s2*xr(459+9*(intsn-1)+jn)* &
                                              antnj(ii, jj)
                end do
            end do
!
        else
!
!------- SUPERNOEUD
!
!---------- PARTIE ROTATION SEULEMENT
!
            do jj = 1, 3
                do ii = 1, 3
                    matn(ii, nb1*6+jj) = -epais*ksi3s2 &
                                         *xr(459+9*(intsn-1)+jn)*antnj( &
                                         ii, jj)
!
                end do
            end do
!
        end if
!
    end do
!
!FIN
!
end subroutine
