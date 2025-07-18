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
subroutine ejcine(ndim, axi, nno1, nno2, vff1, &
                  vff2, wref, dffr2, geom, wg, &
                  kpg, ipg, idf2, rot, b)
! person_in_charge: jerome.laverne at edf.fr
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/r8inir.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "blas/ddot.h"
    aster_logical :: axi
    integer(kind=8) :: ndim, nno1, nno2, kpg, ipg, idf2
    real(kind=8) :: wref, vff1(nno1), vff2(nno2), geom(ndim, nno2)
    real(kind=8) :: rot(ndim, ndim)
    real(kind=8) :: dffr2(ndim-1, nno2), wg, b(2*ndim-1, ndim+1, 2*nno1+nno2)
!
!-----------------------------------------------------------------------
!  MATRICE CINEMATIQUE POUR LES ELEMENTS DE JOINT HM (EN UN PG DONNE)
!-----------------------------------------------------------------------
! IN  NDIM   DIMENSION DE L'ESPACE
! IN  AXI    .TRUE. SI AXISYMETRIQUE
! IN  NNO1   NB DE NOEUDS DE LA FACE POUR LES DEPLACEMENTS
! IN  NNO2   NB DE NOEUDS DE LA FACE POUR LES PRESSIONS P ET LA GEOM X
! IN  VFF1   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR U
! IN  VFF2   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR P ET X
! IN  WREF   POIDS DE REFERENCE DU POINT DE GAUSS
! IN  DFFR2  DERIVEE DES FONCTIONS DE FORME DE REFERENCE DE P ET X EN G
! IN  GEOM   COORDONNEES DES NOEUDS (X)
! OUT WG     POIDS REEL DU POINT DE GAUSS (AVEC DISTORSION)
! OUT ROT    MATRICE DE ROTATION DU LOCAL AU GLOBAL
! OUT B      MATRICE DE PASSAGE UNODAL -> SAUT DE U LOCAL
!-----------------------------------------------------------------------
    integer(kind=8) :: n, i, j
    real(kind=8) :: cova(3, 3), metr(2, 2), cour, jac, cosa, sina, noa1
    real(kind=8) :: ray
    real(kind=8) :: geoloc(ndim, nno2), geotan(ndim-1, nno2)
    real(kind=8) :: dfdis(nno2, ndim-1), wg2
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
    if (ndim .eq. 3) then
!
        call subaco(nno2, dffr2, geom, cova)
        call sumetr(cova, metr, jac)
        wg = wref*jac
!
!       MATRICE DE ROTATION
        noa1 = sqrt(cova(1, 1)**2+cova(2, 1)**2+cova(3, 1)**2)
        rot(1, 1) = cova(1, 3)
        rot(1, 2) = cova(2, 3)
        rot(1, 3) = cova(3, 3)
        rot(2, 1) = cova(1, 1)/noa1
        rot(2, 2) = cova(2, 1)/noa1
        rot(2, 3) = cova(3, 1)/noa1
        rot(3, 1) = rot(1, 2)*rot(2, 3)-rot(1, 3)*rot(2, 2)
        rot(3, 2) = rot(1, 3)*rot(2, 1)-rot(1, 1)*rot(2, 3)
        rot(3, 3) = rot(1, 1)*rot(2, 2)-rot(1, 2)*rot(2, 1)
!
!       CALCUL DE LA GEOMETRIE DANS LE REPERE LOCAL
        call r8inir(ndim*nno2, 0.d0, geoloc, 1)
!
        do n = 1, nno2
            do i = 1, ndim
                do j = 1, ndim
                    geoloc(i, n) = geoloc(i, n)+rot(i, j)*geom(j, n)
                end do
            end do
        end do
!
        do n = 1, nno2
            do i = 2, ndim
                geotan(i-1, n) = geoloc(i, n)
            end do
        end do
!
!       CALCUL DES DERIVEE DES FF DANS LE PLAN TANGENTIEL
        call dfdm2d(nno2, kpg, ipg, idf2, geotan, &
                    wg2, dfdis(1, 1), dfdis(1, 2))
!
    else if (ndim .eq. 2) then
!
!       CALCUL DES DERIVEE DES FF DANS LE PLAN TANGENTIEL
        call dfdm1d(nno2, wref, dffr2, geom, dfdis(1, 1), &
                    cour, wg, cosa, sina)
!
!       CALCUL DE LA DISTANCE A L'AXE EN AXI, R=RAYON DU PG COURANT
        if (axi) then
            b_n = to_blas_int(nno2)
            b_incx = to_blas_int(2)
            b_incy = to_blas_int(1)
            ray = ddot(b_n, geom, b_incx, vff2, b_incy)
            wg = ray*wg
        end if
!
!       MATRICE DE ROTATION
        rot(1, 1) = -cosa
        rot(1, 2) = -sina
        rot(2, 1) = sina
        rot(2, 2) = -cosa
!
    end if
!
!     CONSTRUCTION DE LA MATRICE B
    call r8inir((2*ndim-1)*(ndim+1)*(2*nno1+nno2), 0.d0, b, 1)
!
    do i = 1, ndim
        do j = 1, ndim
!
            do n = 1, nno1
                b(i, j, n) = -rot(i, j)*vff1(n)
                b(i, j, n+nno1) = rot(i, j)*vff1(n)
            end do
!
        end do
    end do
!
    do i = 1, ndim-1
        do n = 1, nno2
            b(ndim+i, ndim+1, 2*nno1+n) = dfdis(n, i)
        end do
    end do
!
end subroutine
