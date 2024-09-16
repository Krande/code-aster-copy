! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine gareac(xdm, xdp, dgamma)
    implicit none
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/angvx.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/promat.h"
#include "asterfort/provec.h"
#include "asterfort/trigom.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
    real(kind=8) :: xdm(3), xdp(3), dgamma
!
    real(kind=8) :: xdmnor(3), xdpnor(3), dd(3)
    real(kind=8) :: normm, normp
    real(kind=8) :: ytr(3), ztr(3)
    real(kind=8) :: ptr2gl(3, 3), pstrx, pstry, ptrbtr(3, 3)
    real(kind=8) :: anglm(3), pglm(3, 3), vtemp(3), pslty, psltz, plo2tr(3, 3)
    real(kind=8) :: ytemp(3), ylocp(3)
    real(kind=8) :: anglp(3), pglp(3, 3)
    real(kind=8) :: pscal, norm
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
! --- SI L'ORIENTATION N'A PAS CHANGE ==> DGAMMA = 0
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, xdm, b_incx, xdmnor, b_incy)
    call normev(xdmnor, normm)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, xdp, b_incx, xdpnor, b_incy)
    call normev(xdpnor, normp)
    dd = xdpnor-xdmnor
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    if (abs(ddot(b_n, dd, b_incx, dd, b_incy)) .lt. r8prem()) then
        dgamma = 0.d0
        goto 9999
    end if
!
! --- CONSTRUCTION DE LA MATRICE DE PASSAGE DU REPERE GLOBAL
!     VERS (XDMNOR,YTR,ZTR)
!     YTR EST NORMAL A XDMNOR DANS LE PLAN FORME PAR (XDMNOR,XDPNOR)
!     ZTR EST NORMAL AU PLAN FORME PAR (XDMNOR,XDPNOR)
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    pscal = ddot(b_n, xdpnor, b_incx, xdmnor, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, xdpnor, b_incx, ytr, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -pscal, xdmnor, b_incx, ytr, &
               b_incy)
    call normev(ytr, norm)
    if (norm .le. r8miem()) then
        dgamma = 0.d0
        goto 9999
    end if
    call provec(xdmnor, ytr, ztr)
!
    ptr2gl(1, 1) = xdmnor(1)
    ptr2gl(2, 1) = xdmnor(2)
    ptr2gl(3, 1) = xdmnor(3)
    ptr2gl(1, 2) = ytr(1)
    ptr2gl(2, 2) = ytr(2)
    ptr2gl(3, 2) = ytr(3)
    ptr2gl(1, 3) = ztr(1)
    ptr2gl(2, 3) = ztr(2)
    ptr2gl(3, 3) = ztr(3)
!
! --- CONSTRUCTION DE LA MATRICE DE PASSAGE DE (XDMNOR,YTR,ZTR)
!     VERS (XDPNOR,YDPNOR,ZTR)
!     IL S'AGIT DE LA ROTATION AUTOUR DE ZTR
!     YDPNOR EST NORMAL A XDPNOR DANS LE PLAN NORMAL A ZTR
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    pstrx = ddot(b_n, xdmnor, b_incx, xdpnor, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    pstry = ddot(b_n, ytr, b_incx, xdpnor, b_incy)
!
    ptrbtr(1, 1) = pstrx
    ptrbtr(1, 2) = -pstry
    ptrbtr(1, 3) = 0.d0
    ptrbtr(2, 1) = pstry
    ptrbtr(2, 2) = pstrx
    ptrbtr(2, 3) = 0.d0
    ptrbtr(3, 1) = 0.d0
    ptrbtr(3, 2) = 0.d0
    ptrbtr(3, 3) = 1.d0
!
! --- CONSTRUCTION DE LA MATRICE DE PASSAGE DE (XDMNOR,YTR,ZTR)
!     VERS (XDMNOR,YLOCM,ZLOCM)
!     IL S'AGIT DE LA ROTATION AUTOUR DE XDMNOR
!     (YLOCM,ZLOCM) CONSTITUE LE REPERE LOCAL A T- POUR GAMMA=0
!
    call angvx(xdm, anglm(1), anglm(2))
    anglm(3) = 0.d0
    call matrot(anglm, pglm)
    vtemp(1) = pglm(2, 1)
    vtemp(2) = pglm(2, 2)
    vtemp(3) = pglm(2, 3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    pslty = ddot(b_n, ytr, b_incx, vtemp, b_incy)
    vtemp(1) = pglm(3, 1)
    vtemp(2) = pglm(3, 2)
    vtemp(3) = pglm(3, 3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    psltz = ddot(b_n, ytr, b_incx, vtemp, b_incy)
!
    plo2tr(1, 1) = 1.d0
    plo2tr(1, 2) = 0.d0
    plo2tr(1, 3) = 0.d0
    plo2tr(2, 1) = 0.d0
    plo2tr(2, 2) = pslty
    plo2tr(2, 3) = psltz
    plo2tr(3, 1) = 0.d0
    plo2tr(3, 2) = -psltz
    plo2tr(3, 3) = pslty
!
! --- CALCUL DE YLOC A T+
!
!     PTR2GL * PTRBTR * PLO2TR *  (0 1 0)T
    vtemp(1) = plo2tr(1, 2)
    vtemp(2) = plo2tr(2, 2)
    vtemp(3) = plo2tr(3, 2)
    call promat(ptrbtr, 3, 3, 3, vtemp, &
                3, 3, 1, ytemp)
    call promat(ptr2gl, 3, 3, 3, ytemp, &
                3, 3, 1, ylocp)
!
! --- CALCUL DE LA BASE LOCALE A T+ EN CONSIDERANT GAMMA NUL
    call angvx(xdp, anglp(1), anglp(2))
    anglp(3) = 0.d0
    call matrot(anglp, pglp)
!     LES LIGNES DE PGLP SONT LES 3 VECTEURS DE LA BASE LOCAL
!
    vtemp(1) = pglp(2, 1)
    vtemp(2) = pglp(2, 2)
    vtemp(3) = pglp(2, 3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    pscal = ddot(b_n, ylocp, b_incx, vtemp, b_incy)
!
    dgamma = trigom('ACOS', pscal)
    vtemp(1) = pglp(3, 1)
    vtemp(2) = pglp(3, 2)
    vtemp(3) = pglp(3, 3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    if (ddot(b_n, ylocp, b_incx, vtemp, b_incy) .lt. 0.d0) then
        dgamma = -dgamma
    end if
!
9999 continue
!
end subroutine
