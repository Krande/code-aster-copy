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
subroutine ejcine_hm(ndim, axi, nno1, nno2, vff1, &
                     vff2, wref, dffr2, geom, ang, &
                     wg, b)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/matrot.h"
#include "asterfort/r8inir.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "blas/ddot.h"
    aster_logical :: axi
    integer(kind=8) :: ndim, nno1, nno2
    real(kind=8) :: wref, vff1(nno1), vff2(nno2), geom(ndim, nno2), ang(*)
    real(kind=8) :: dffr2(ndim-1, nno2), wg, b(3, 3, 2*nno1)
!-----------------------------------------------------------------------
!  MATRICE CINEMATIQUE POUR LES ELEMENTS D'INTERFACE (EN UN PG DONNE)
!-----------------------------------------------------------------------
! IN  NDIM   DIMENSION DE L'ESPACE
! IN  AXI    .TRUE. SI AXISYMETRIQUE
! IN  NNO1   NB DE NOEUDS DE LA FACE POUR LES DEPLACEMENTS
! IN  NNO2   NB DE NOEUDS DE LA FACE POUR LES LAGRANGES L ET LA GEOM X
! IN  VFF1   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR U
! IN  VFF2   VALEUR DES FONCTIONS DE FORME (DE LA FACE) POUR L ET X
! IN  WREF   POIDS DE REFERENCE DU POINT DE GAUSS
! IN  DFFR2  DERIVEE DES FONCTIONS DE FORME DE REFERENCE DE L ET X EN G
! IN  GEOM   COORDONNEES DES NOEUDS (X)
! IN  ANG    ANGLES NAUTIQUES NODAUX (FAMILLE X)
! OUT WG     POIDS REEL DU POINT DE GAUSS (AVEC DISTORSION)
! OUT B      MATRICE DE PASSAGE UNODAL -> SAUT DE U LOCAL
!-----------------------------------------------------------------------
    integer(kind=8) :: n, i, j, nang
    real(kind=8) :: cova(3, 3), metr(2, 2), dfdx(9), cour, jac, cosa, sina
    real(kind=8) :: angloc(3), rot(3, 3), r, rmax
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
!    CALCUL DU JACOBIEN
!
    if (ndim .eq. 3) then
        call subaco(nno2, dffr2, geom, cova)
        call sumetr(cova, metr, jac)
        wg = wref*jac
    else if (ndim .eq. 2) then
        call dfdm1d(nno2, wref, dffr2, geom, dfdx, &
                    cour, wg, cosa, sina)
    end if
!
    if (axi) then
        b_n = to_blas_int(nno2)
        b_incx = to_blas_int(2)
        b_incy = to_blas_int(1)
        r = ddot(b_n, geom, b_incx, vff2, b_incy)
! ----------------------------------------------------------------------
! POUR LES ELEMENTS AVEC COUPLAGE HM, DANS LE CAS OU R EGAL 0, ON A UN
! JACOBIEN NUL EN UN PG. ON PRENDS LE MAX DU RAYON MULTIPLIE PAR 1.E-3
! ----------------------------------------------------------------------
        if (r .eq. 0.d0) then
            rmax = geom(1, 1)
            do n = 2, nno2
                rmax = max(geom(1, n), rmax)
            end do
            wg = wg*1.d-03*rmax
        else
            wg = r*wg
        end if
    end if
!
!    CALCUL DES ANGLES NAUTIQUES AU POINT D'INTEGRATION
!
    if (ndim .eq. 2) nang = 1
    if (ndim .eq. 3) nang = 3
    call r8inir(3, 0.d0, angloc, 1)
    do i = 1, nang
        b_n = to_blas_int(nno2)
        b_incx = to_blas_int(nang)
        b_incy = to_blas_int(1)
        angloc(i) = ddot(b_n, ang(i), b_incx, vff2, b_incy)
    end do
!
!    CALCUL DE LA MATRICE DE ROTATION GLOBAL -> LOCAL
!
    call matrot(angloc, rot)
!
!    CONSTRUCTION DE LA MATRICE B
!
    do i = 1, ndim
        do j = 1, ndim
            do n = 1, nno1
                b(i, j, n) = -rot(i, j)*vff1(n)
                b(i, j, n+nno1) = rot(i, j)*vff1(n)
            end do
        end do
    end do
!
end subroutine
