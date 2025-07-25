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
subroutine symema(geomi, perp, pt)
    implicit none
!
!     BUT : SYMETRIE D'UN MAILLAGE PAR RAPPORT A UN PLAN
!
!     IN :
!            GEOMI  : CHAM_NO(GEOM_R) : CHAMP DE GEOMETRIE A SYMETRISER
!            PERP   : AXE PERPENDICULAIRE AU PLAN
!            PT     : UN POINT DU PLAN
!     OUT:
!            GEOMI  : CHAM_NO(GEOM_R) : CHAMP DE GEOMETRIE ACTUALISE
!
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
    character(len=24) :: coorjv
    character(len=19) :: geomi
    real(kind=8) :: norm, prec, xd, pti(3), pt(3), perp(3), dist
    integer(kind=8) :: i, iadcoo, n1
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call matfpe(-1)
!
    call jemarq()
!     RECUPERATION DE L'ADRESSE DES COORDONNEES ET DU NOMBRE DE POINTS
    coorjv = geomi(1:19)//'.VALE'
    call jeveuo(coorjv, 'E', iadcoo)
    call jelira(coorjv, 'LONMAX', n1)
    iadcoo = iadcoo-1
    n1 = n1/3
!
!     NORMALISATION DE PERP
    prec = 1.d-14
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    norm = dnrm2(b_n, perp, b_incx)
    if (norm .lt. prec) then
        call utmess('F', 'ALGORITH10_87')
    end if
    perp(1) = perp(1)/norm
    perp(2) = perp(2)/norm
    perp(3) = perp(3)/norm
!     LE PLAN PASSE PAR "PT"
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    xd = -ddot(b_n, perp, b_incx, pt, b_incy)
!
!     BOUCLE SUR TOUS LES POINTS
    do i = 1, n1
        pti(1) = zr(iadcoo+3*(i-1)+1)
        pti(2) = zr(iadcoo+3*(i-1)+2)
        pti(3) = zr(iadcoo+3*(i-1)+3)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        dist = ddot(b_n, perp, b_incx, pti, b_incy)+xd
        zr(iadcoo+3*(i-1)+1) = -2.0d0*dist*perp(1)+pti(1)
        zr(iadcoo+3*(i-1)+2) = -2.0d0*dist*perp(2)+pti(2)
        zr(iadcoo+3*(i-1)+3) = -2.0d0*dist*perp(3)+pti(3)
    end do
!
    call jedema()
    call matfpe(1)
!
end subroutine
