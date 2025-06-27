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
subroutine porea2(nno, nc, geom, gamma, pgl, &
                  xl, dispParaNameZ_)
    implicit none
#include "jeveux.h"
#include "asterfort/angvx.h"
#include "asterfort/assert.h"
#include "asterfort/matrot.h"
#include "asterfort/tecach.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: nno, nc
    real(kind=8) :: geom(3, nno), gamma
!
    real(kind=8) :: pgl(3, 3), xl
    character(len=*), optional, intent(in) :: dispParaNameZ_
!
!
!     ------------------------------------------------------------------
!     CALCUL DE LA MATRICE DE PASSAGE GLOBALE/LOCALE EN TENANT COMPTE
!     DE LA GEOMETRIE REACTUALISEE AINSI QUE LA LONGUEUR DE LA POUTRE
!     POUR L'OPTION FORC_NODA
!     ------------------------------------------------------------------
! IN  NNO    : NOMBRE DE NOEUDS
! IN  NC     : NOMBRE DE COMPOSANTE DU CHAMP DE DEPLACEMENTS
! IN  GEOM   : COORDONNEES DES NOEUDS
! IN  GAMMA  : ANGLE DE VRILLE AU TEMPS -
! OUT PGL    : MATRICE DE PASSAGE GLOBAL/LOCAL
!     ------------------------------------------------------------------
!
!     VARIABLES LOCALES
    integer(kind=8) :: i, ideplm, ideplp, iret
    real(kind=8) :: utg(14), xug(6), xd(3), xl2, alfa1, beta1
    real(kind=8) :: tet1, tet2, gamma1, ang1(3)
    character(len=16) :: dispParaName
    blas_int :: b_incx, b_incy, b_n
!
    ASSERT(nno .eq. 2)
    dispParaName = "PDEPLMR"
    if (present(dispParaNameZ_)) then
        dispParaName = dispParaNameZ_
    end if
!
    call tecach('ONO', dispParaName, 'L', iret, iad=ideplm)
    if (iret .ne. 0) then
        utg = 0.d0
    else
        do i = 1, 2*nc
            utg(i) = zr(ideplm-1+i)
        end do
    end if
!
    if (dispParaName .ne. "PDEPLAR") then
        call tecach('ONO', 'PDEPLPR', 'L', iret, iad=ideplp)
        if (iret .eq. 0) then
            do i = 1, 2*nc
                utg(i) = utg(i)+zr(ideplp-1+i)
            end do
        end if
    end if
!
    do i = 1, 3
        xug(i) = utg(i)+geom(i, 1)
        xug(i+3) = utg(i+nc)+geom(i, 2)
        xd(i) = xug(i+3)-xug(i)
    end do
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    xl2 = ddot(b_n, xd, b_incx, xd, b_incy)
    xl = sqrt(xl2)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    tet1 = ddot(b_n, utg(4), b_incx, xd, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    tet2 = ddot(b_n, utg(nc+4), b_incx, xd, b_incy)
    tet1 = tet1/xl
    tet2 = tet2/xl
    call angvx(xd, alfa1, beta1)
!
    gamma1 = gamma+(tet1+tet2)/2.d0
    ang1(1) = alfa1
    ang1(2) = beta1
    ang1(3) = gamma1
!
    call matrot(ang1, pgl)
end subroutine
