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
subroutine chgref(geomi, x, y, bidim)
    implicit none
!
!     BUT : CHANGEMENT DE REPERE D'UN MAILLAGE
!
!     IN :
!            GEOMI  : CHAM_NO(GEOM_R) : CHAMP DE GEOMETRIE A TOURNER
!            X      : PREMIER VECTEUR DE LA BASE
!            Y      : DEUXIEME VECTEUR DE LA BASE
!            BIDIM  : BOOLEEN VRAI SI GEOMETRIE 2D
!     OUT:
!            GEOMI  : CHAM_NO(GEOM_R) : CHAMP DE GEOMETRIE ACTUALISE
!
!
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
    integer(kind=8) :: n1, i, iadcoo
    aster_logical :: bidim
    character(len=19) :: geomi
    character(len=24) :: coorjv
    real(kind=8) :: x(3), y(3), z(3), p(3), prec, r1, r2
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
!
    call matfpe(-1)
!
    call jemarq()
    coorjv = geomi(1:19)//'.VALE'
    call jeveuo(coorjv, 'E', iadcoo)
    call jelira(coorjv, 'LONMAX', n1)
    prec = 1.d-14
    n1 = n1/3
    iadcoo = iadcoo-1
!     -- ON TRAITE LE CAS 2D SEPAREMENT POUR OPTIMISER :
    if (bidim) then
        b_n = to_blas_int(2)
        b_incx = to_blas_int(1)
        r1 = dnrm2(b_n, x, b_incx)
        if (r1 .gt. prec) then
            x(1) = x(1)/r1
            x(2) = x(2)/r1
            x(3) = x(3)/r1
            y(1) = (-1.d0)*x(2)
            y(2) = x(1)
            do i = 1, n1
                p(1) = zr(iadcoo+3*(i-1)+1)
                p(2) = zr(iadcoo+3*(i-1)+2)
                b_n = to_blas_int(2)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iadcoo+3*(i-1)+1) = ddot(b_n, x, b_incx, p, b_incy)
                b_n = to_blas_int(2)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iadcoo+3*(i-1)+2) = ddot(b_n, y, b_incx, p, b_incy)
            end do
        else
            call utmess('F', 'ALGORITH_96')
        end if
    else
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        r1 = dnrm2(b_n, x, b_incx)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        r2 = dnrm2(b_n, y, b_incx)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if ((ddot(b_n, x, b_incx, y, b_incy) .lt. prec) .and. ((r1*r2) .gt. 0.d0)) then
            x(1) = x(1)/r1
            x(2) = x(2)/r1
            x(3) = x(3)/r1
            y(1) = y(1)/r2
            y(2) = y(2)/r2
            y(3) = y(3)/r2
            call provec(x, y, z)
            do i = 1, n1
                p(1) = zr(iadcoo+3*(i-1)+1)
                p(2) = zr(iadcoo+3*(i-1)+2)
                p(3) = zr(iadcoo+3*(i-1)+3)
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iadcoo+3*(i-1)+1) = ddot(b_n, x, b_incx, p, b_incy)
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iadcoo+3*(i-1)+2) = ddot(b_n, y, b_incx, p, b_incy)
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iadcoo+3*(i-1)+3) = ddot(b_n, z, b_incx, p, b_incy)
            end do
        else
            call utmess('F', 'ALGORITH_97')
        end if
    end if
    call jedema()
    call matfpe(1)
!
end subroutine
