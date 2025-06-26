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
subroutine arcseg34(nbno, coor, abscur)
!
    implicit none
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/assert.h"
#include "blas/ddot.h"
#include "asterfort/mgauss.h"
#include "asterfort/utmess.h"
#include "asterfort/provec.h"
!
    integer(kind=8), intent(in) :: nbno
    real(kind=8), intent(in) :: coor(3, nbno)
    real(kind=8), intent(out) :: abscur(nbno)
! ......................................................................
!  but:  calcul de l'abscisse curviligne de 3 ou 4 points situés sur
!        un arc de cercle
!
!  arguments :
!     nbno : nombre de noeuds de l'arc (SEG3 ou SEG4)
!     coor : coordonnees des noeuds de l'arc
!     abscur : abscisse curviligne des noeuds sur l'arc
!
!  remarques importantes :
!     1) l'arc est suppose etre proche d'un arc de cercle.
!     2) l'ordre des noeuds est celui des SEG3/4 :
!              1 - 3 - 2       (SEG3)
!        ou    1 - 3 - 4 - 2   (SEG4)
!     3) => abscur(1)=0, abscur(2)=longueur totale de l'arc
! ......................................................................
!
    real(kind=8) :: a(3), b(3), c(3), ab(3), bc(3), ce(3), mab(3), mbc(3)
    real(kind=8) :: n(3), mat(3, 3), r, r2, ra(3), rk(3), x(3), det
    real(kind=8) :: sintheta, costheta, theta, valr(6)
    integer(kind=8) :: k, iret
    blas_int :: b_incx, b_incy, b_n
!   ----------------------------------------------------------------------------
    ASSERT(nbno .ge. 1 .and. nbno .le. 4)
!
!   Soit A, B et C 3 noeuds de l'arc (A et C sont les extremites)
    a(:) = coor(:, 1)
    c(:) = coor(:, 2)
    b(:) = coor(:, 3)
    ab(:) = b(:)-a(:)
    bc(:) = c(:)-b(:)
    call provec(ab, bc, n)
!
!   -- 1. Cas des SEG2 (forcement droit) :
!   ---------------------------------------
    if (nbno .eq. 2) then
        abscur(1) = 0.d0
        x(:) = coor(:, 2)-coor(:, 1)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        abscur(2) = abscur(1)+sqrt(ddot(b_n, x, b_incx, x, b_incy))
        goto 999
    end if
!
!   -- 1. Cas des points presque alignes :
!   ---------------------------------------
!   -- les points sont presqu'alignes si l'angle (ab,bc) < 1 degre :
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    sintheta = sqrt( &
               ddot(b_n, n, b_incx, n, b_incy))/(sqrt(ddot(b_n, ab, b_incx, ab, b_incy))*sqrt(ddo&
               &t(b_n, bc, b_incx, bc, b_incy)) &
               )
    theta = abs(asin(sintheta))
    if ((theta*r8rddg()) .lt. 1.) then
!       -- calcul des abscisses curvilignes :
        abscur(1) = 0.d0
        x(:) = coor(:, 3)-coor(:, 1)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        abscur(3) = abscur(1)+sqrt(ddot(b_n, x, b_incx, x, b_incy))
        if (nbno .eq. 3) then
            x(:) = coor(:, 2)-coor(:, 3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            abscur(2) = abscur(3)+sqrt(ddot(b_n, x, b_incx, x, b_incy))
        else
            x(:) = coor(:, 4)-coor(:, 3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            abscur(4) = abscur(3)+sqrt(ddot(b_n, x, b_incx, x, b_incy))
            x(:) = coor(:, 2)-coor(:, 4)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            abscur(2) = abscur(4)+sqrt(ddot(b_n, x, b_incx, x, b_incy))
        end if
        goto 999
    end if
!
!   -- 2. Cas des points formant un arc de cercle :
!   ------------------------------------------------
!
!   -- 2.1 Determination du centre du cercle avec les 3 premiers noeuds (CE) :
!   Le centre CE est l'intersection des 3 plans :
!     * P1 : le plan mediateur de AB
!     * P2 : le plan mediateur de BC
!     * P3 : le plan contenant A, B et C
    mab(:) = (a(:)+b(:))/2.d0
    mbc(:) = (b(:)+c(:))/2.d0
!
!   -- mat est la matrice (3,3) contenant les coefs. des equations des 3 plans :
!   -- ce  est le vecteur (3) contenant les termes constants des 3 plans :
!
!   -- plan P1 :
    mat(1, :) = ab(:)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ce(1) = ddot(b_n, mab, b_incx, ab, b_incy)
!
!   -- plan P2 :
    mat(2, :) = bc(:)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ce(2) = ddot(b_n, mbc, b_incx, bc, b_incy)
!
!   -- plan P3 :
    mat(3, :) = n(:)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ce(3) = ddot(b_n, n, b_incx, a, b_incy)
!
!   -- resolution de mat*ce=ce pour trouver CE :
    call mgauss('NFSP', mat, ce, 3, 3, &
                1, det, iret)
    ASSERT(iret .eq. 0)
!
!   -- on verifie que les points sont bien sur le cercle trouve (a 1% pres):
    r = sqrt((ce(1)-a(1))**2+(ce(2)-a(2))**2+(ce(3)-a(3))**2)
    ASSERT(r .gt. 0.d0)
    do k = 2, nbno
        r2 = sqrt((ce(1)-coor(1, k))**2+(ce(2)-coor(2, k))**2+(ce(3)-coor(3, k))**2)
        if (abs((r2-r)/r) .gt. 1.d-2) then
            valr(1:3) = a(:)
            valr(4:6) = c(:)
            call utmess('F', 'MODELISA_4', nr=6, si=nbno, valr=valr)
        end if
    end do
!
!   -- calcul des abscisses curvilignes :
    abscur(1) = 0.d0
    ra(:) = ce(:)-a(:)
    do k = 2, nbno
        rk(:) = ce(:)-coor(:, k)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        costheta = ddot(b_n, ra, b_incx, rk, b_incy)/(r*r)
        theta = acos(costheta)
        abscur(k) = r*theta
    end do
!
999 continue
!
end subroutine
