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
subroutine lcresa(fami, kpg, ksp, typmod, imat, &
                  nmat, materd, materf, rela_comp, nr, &
                  nvi, timed, timef, deps, epsd, &
                  yf, dy, r, iret, yd, &
                  crit)
!
! aslint: disable=W1306,W1504
    implicit none
!     CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = R(DY)
!     IN  FAMI   :  FAMILLE DU POINT DE GAUSS
!         KPG    :  POINT DE GAUSS
!         KSP    :  SOUS-POINT DE GAUSS
!         LOI    :  MODELE DE COMPORTEMENT
!         TYPMOD    :  TYPE DE MODELISATION
!         IMAT   :  NOM DU MATERIAU
!         NMAT   :  DIMENSION MATER
!         MATERD :  COEFFICIENTS MATERIAU A T
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!         COMP   :  COMPORTEMENT
!         NR     :  DIMENSION DU SYSTEME R(NR) NR=6+NVI
!         NVI    :  NOMBRE DE VARIABLE INTERNES
!         TIMED  :  INSTANT  T
!         TIMEF  :  INSTANT  T+DT
!         DEPS   :  INCREMENT DE DEFORMATION
!         EPSD   :  DEFORMATION A T
!         YF     :  VARIABLES A T + DT
!         DY     :  SOLUTION
!     OUT R      :  SYSTEME NL A T + DT
!     ----------------------------------------------------------------
!
#include "asterfort/calsig.h"
#include "asterfort/lcdvin.h"
#include "asterfort/lcopil.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: imat, nmat, nr, nvi, kpg, ksp, iret, itens, ndt, ndi
    real(kind=8) :: deps(6), epsd(6), r(nr), yf(nr), dy(nr), x, theta, yd(*)
    real(kind=8) :: crit(*)
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), timed, timef, evi(6)
    real(kind=8) :: dkooh(6, 6), fkooh(6, 6), h1sigf(6), sigi(6), vini(nvi)
    real(kind=8) :: dtime, dvin(nvi), epsef(6), sigf(6), smx(6), dsig(6)
    character(len=8) :: typmod
    character(len=*) :: fami
    character(len=3) :: matcst
    character(len=16) :: rela_comp
    blas_int :: b_incx, b_incy, b_n
    common/tdim/ndt, ndi
!----------------------------------------------------------------
!
    iret = 0
    dtime = timef-timed
    theta = crit(4)
!
!     INVERSE DE L'OPERATEUR D'ELASTICITE DE HOOKE
    if (materf(nmat, 1) .eq. 0) then
        call lcopil('ISOTROPE', typmod, materd, dkooh)
        call lcopil('ISOTROPE', typmod, materf, fkooh)
    else if (materf(nmat, 1) .eq. 1) then
        call lcopil('ORTHOTRO', typmod, materd, dkooh)
        call lcopil('ORTHOTRO', typmod, materf, fkooh)
    end if
!
    x = dtime
    matcst = 'OUI'
!
!     CALCUL DE DSIG AVEC THETA
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, dy, b_incx, dsig, b_incy)
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, yd, b_incx, smx, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, theta, dsig, b_incx, smx, &
               b_incy)
!
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, yd(ndt+1), b_incx, vini, b_incy)
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, theta, dy, b_incx, vini, &
               b_incy)
!
!     CALCUL DES DERIVEES DES VARIABLES INTERNES AU POINT T+THETA*DT
    call lcdvin(fami, kpg, ksp, rela_comp, typmod, &
                imat, matcst, nvi, nmat, vini, &
                materf(1, 2), x, dtime, smx, dvin, &
                iret)
!
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    call dscal(b_n, dtime, dvin, b_incx)
!
    do itens = 1, nvi
        vini(itens) = yd(ndt+itens)+dvin(itens)
    end do
!
    do itens = 1, 6
        evi(itens) = vini(itens)
    end do
!
!     CALCUL DES CONTRAINTES AU POINT T+DT
    call calsig(fami, kpg, ksp, evi, typmod, &
                rela_comp, vini, x, dtime, epsd, &
                deps, nmat, materf(1, 1), sigi)
!
!     CALCUL DES RESIDUS AU POINT T+DT
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, yf, b_incx, sigf, b_incy)
    epsef(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigi(1:ndt))
    h1sigf(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
    r(1:ndt) = epsef(1:ndt)-h1sigf(1:ndt)
!
!     CALCUL DES RESIDUS AU POINT T+DT
!
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, dvin, b_incx, r(ndt+1), b_incy)
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, dy(ndt+1), b_incx, r(ndt+1), &
               b_incy)
!
end subroutine
