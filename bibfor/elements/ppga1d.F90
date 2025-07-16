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
subroutine ppga1d(ndim, nno, npg, poids, vff, &
                  dff, geom, pg)
    implicit none
#include "jeveux.h"
    integer(kind=8) :: ndim, nno, npg
    real(kind=8) :: poids(npg), vff(nno, npg), dff(nno, npg), geom(ndim, nno)
    real(kind=8) :: pg(ndim+1, npg)
! ----------------------------------------------------------------------
!  POSITION ET POIDS DES POINTS DE GAUSS POUR ELEMENTS 1D
!  *           *                   **                  **
! ----------------------------------------------------------------------
! IN  NDIM    DIMENSION DE L'ESPACE
! IN  NNO     NOMBRE DE NOEUDS
! IN  NPG     NOMBRE DE POINTS DE GAUSS
! IN  POIDS   POIDS DES POINTS DE GAUSS DE L'ELEMENT DE REFERENCE
! IN  VFF     VALEUR DES FONCTIONS DE FORME
! IN  GEOM    COORDONNEES DES NOEUDS
! OUT PG      COORDONNEES DES POINTS DE GAUSS + POIDS
! ----------------------------------------------------------------------
!
    integer(kind=8) :: g, i, j
    real(kind=8) :: dxdk, dydk, dzdk
    real(kind=8) :: jac
!
! ----------------------------------------------------------------------
!
!     1. CALCUL DES COORDONNEES DES POINTS DE GAUSS
!     =============================================
!
    do g = 1, npg
        do i = 1, ndim+1
            pg(i, g) = 0.d0
        end do
    end do
!
    do g = 1, npg
        do i = 1, ndim
            do j = 1, nno
                pg(i, g) = pg(i, g)+geom(i, j)*vff(j, g)
            end do
        end do
    end do
!
    if (nno .eq. 1) goto 999
!
!     2. CALCUL DU POIDS
!     ==================
!
    do g = 1, npg
        dxdk = 0.d0
        dydk = 0.d0
        dzdk = 0.d0
!        COMPOSANTES DU VECTEUR TANGENT
        do j = 1, nno
            dxdk = dxdk+geom(1, j)*dff(j, g)
            dydk = dydk+geom(2, j)*dff(j, g)
            if (ndim .eq. 3) dzdk = dzdk+geom(3, j)*dff(j, g)
        end do
!        JACOBIEN 1D == DERIVEE DE L'ABSCISSE CURVILIGNE
        jac = sqrt(dxdk**2+dydk**2+dzdk**2)
!
        pg(ndim+1, g) = jac*poids(g)
    end do
!
999 continue
!
end subroutine
