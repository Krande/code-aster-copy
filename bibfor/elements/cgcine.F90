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
subroutine cgcine(ndim, nno1, vff1, wref, dffr1, &
                  geom, tang, wg, l, b, &
                  nornor)
!
    implicit none
#include "asterfort/dfdm1b.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nno1, ndim
    real(kind=8) :: wref, vff1(nno1), geom(ndim, nno1), tang(3, 3)
    real(kind=8) :: dffr1(nno1), wg, dfdx(3), b(4, nno1), l(nno1)
!-----------------------------------------------------------------------
!  MATRICE CINEMATIQUE POUR LES ELEMENTS DE CABLE/GAINE (EN UN PG DONNE)
!     ROUTINE INSPIREE DE EICINE
!-----------------------------------------------------------------------
! IN  NDIM   DIMENSION DE L'ESPACE
! IN  NNO1   NB DE NOEUDS POUR LES DEPLACEMENTS U
! IN  VFF1   VALEUR DES FONCTIONS DE FORME POUR U
! IN  WREF   POIDS DE REFERENCE DU POINT DE GAUSS
! IN  DFFR1  DERIVEE DES FONCTIONS DE FORME DE REFERENCE DE L EN G
! IN  GEOM   COORDONNEES DES NOEUDS (X)
! IN  TANG   TANGENTES (FAMILLE X)
! OUT WG     POIDS REEL DU POINT DE GAUSS (AVEC DISTORSION)
! OUT L      VALEURS DES FONCTIONS DE FORME
! OUT B      MATRICE DE PASSAGE UNODAL -> U GAINE TANGENTIEL ET U CABLE
! OUT NORNOR RAYON DE COURBURE
!-----------------------------------------------------------------------
    integer(kind=8) :: n, i
    real(kind=8) :: tanloc(3), norloc(3)
    real(kind=8) :: norm, nornor
!-----------------------------------------------------------------------
!
!    CALCUL DU JACOBIEN
!
    call dfdm1b(nno1, wref, dffr1, geom, dfdx, &
                wg)
!
!
!    CALCUL DES ANGLES NAUTIQUES AU POINT D'INTEGRATION
!
    call r8inir(3, 0.d0, tanloc, 1)
    call r8inir(3, 0.d0, norloc, 1)
    do n = 1, nno1
        do i = 1, 3
            tanloc(i) = tanloc(i)+vff1(n)*tang(i, n)
            norloc(i) = norloc(i)+dfdx(n)*tang(i, n)
        end do
    end do
    norm = sqrt((tanloc(1)**2+tanloc(2)**2+tanloc(3)**2))
    nornor = sqrt((norloc(1)**2+norloc(2)**2+norloc(3)**2))
!
!    CONSTRUCTION DE LA MATRICE B
!
    do n = 1, nno1
        do i = 1, ndim
            b(i, n) = tanloc(i)/norm*dfdx(n)
        end do
        b(ndim+1, n) = dfdx(n)
        l(n) = vff1(n)
    end do
!
end subroutine
