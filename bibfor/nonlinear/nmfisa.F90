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

subroutine nmfisa(axi, geom, kpg, poids, b)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/r8inir.h"
    aster_logical :: axi
    integer(kind=8) :: kpg
    real(kind=8) :: geom(2, 4), poids, b(2, 8)
!
!-----------------------------------------------------------------------
!
! BUT:
!     POUR LE POINT DE GAUSS KPG :
!     CALCUL DU POIDS DU POINT DE GAUSS
!     CALCUL DE LA MATRICE B DONNANT LES SAUT PAR ELEMENTS A PARTIR DES
!     DEPLACEMENTS AUX NOEUDS : SU = B U
!
! REMARQUE :
!
!   LA MATRICE B INCLUE LE CHANGEMENT DE REPERE LOCAL/GLOBAL :
!     U    = DEPLACEMENT DANS LE REPERE GLOBAL
!     ULOC = DEPLACEMENT DANS LE REPERE LOCAL A L'ELEMENT
!
!   B S'ECRIT SOUS LA FORME : B = BTILD RTILD
!   AVEC :
!            SU   = BTILD ULOC
!            ULOC = RTILD U
!
!
! IN  : AXI TRUE EN AXI
! IN  : GEOM,KPG
! OUT : POIDS, B
!
!-----------------------------------------------------------------------
    real(kind=8) :: co, si, c, s, coef(2), aire, rayon
!-----------------------------------------------------------------------
!
    coef(1) = 0.5d0*(1.d0+sqrt(3.d0)/3.d0)
    coef(2) = 0.5d0*(1.d0-sqrt(3.d0)/3.d0)
!
! -- CALCUL DU POIDS DU POINT DE GAUSS
!
    aire = sqrt((geom(1, 2)-geom(1, 1))**2+(geom(2, 2)-geom(2, 1))**2)
!
!     DISTANCE DU PG A L'AXE DE REVOLUTION (I.E ABSCISSE DU PG)
    rayon = geom(1, 1)*coef(kpg)+geom(1, 2)*(1.d0-coef(kpg))
    if (axi) aire = aire*rayon
!
    poids = aire/2
!
! -- CALCUL DE LA MATRICE B
!
    call r8inir(16, 0.d0, b, 1)
!
    co = (geom(2, 2)-geom(2, 1))
    si = -(geom(1, 2)-geom(1, 1))
!
    c = co/sqrt(co*co+si*si)
    s = si/sqrt(co*co+si*si)
!
! SAISIE DE LA MATRICE B : APPLICATION LINEAIRE DONNANT LE SAUT DE
! DEPLACEMENT DANS L'ELEMENT (SU_N,SU_T) A PARTIR DES DEPLACEMENTS
! AUX NOEUDS :
!
    b(1, 1) = c*coef(kpg)
    b(1, 2) = s*coef(kpg)
    b(1, 3) = c*(1.d0-coef(kpg))
    b(1, 4) = s*(1.d0-coef(kpg))
    b(1, 5) = -c*(1.d0-coef(kpg))
    b(1, 6) = -s*(1.d0-coef(kpg))
    b(1, 7) = -c*coef(kpg)
    b(1, 8) = -s*coef(kpg)
!
    b(2, 1) = -s*coef(kpg)
    b(2, 2) = c*coef(kpg)
    b(2, 3) = -s*(1.d0-coef(kpg))
    b(2, 4) = c*(1.d0-coef(kpg))
    b(2, 5) = s*(1.d0-coef(kpg))
    b(2, 6) = -c*(1.d0-coef(kpg))
    b(2, 7) = s*coef(kpg)
    b(2, 8) = -c*coef(kpg)
!
end subroutine
