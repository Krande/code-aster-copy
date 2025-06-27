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
subroutine radial(nbsig, sigm, sigp, indm, indp, &
                  icine, xm, xp, normdn)
!
!     BUT:
!       CALCUL DE L'INDICATEUR LOCAL DE RADIALITE
!       NORME DE LA DIFFERENCE DES NORMALES
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NBSIG    : CDIMENSION DES VECTEURS CONTRAINTES
! IN   SIGM     : CONTRAINTES INSTANT + AVEC RAC2
! IN   SIGP     : CONTRAINTES INSTANT - AVEC RAC2
! IN   INDM     : INDICATEUR PLASTICITE INSTANT -
! IN   INDP     : INDICATEUR PLASTICITE INSTANT +
! IN   ICINE    : INDICATEUR ECROUISSAGE CINAMATIQUE
! IN   XM       : VARIABLES INTERNES X INSTANT - AVEC RAC2
! IN   XP       : VARIABLES INTERNES X INSTANT + AVEC RAC2
!
!      SORTIE :
!-------------
! OUT  NORMDN   : INDICATEUR DE PERTE DE RADIALITE=NORME DE DN
!
! ......................................................................
    implicit none
#include "asterc/r8prem.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    integer(kind=8) :: nbsig, icine, i
    real(kind=8) :: n1(6), n2(6), xm(6), xp(6), indm, indp
    real(kind=8) :: tensm(6), tensp(6), normdn, sigm(nbsig), sigp(nbsig)
    real(kind=8) :: zernor, devm(6), devp(6), smeq, speq
    real(kind=8) :: dn1n2(6), ndn, preci, trt1, trt2
    blas_int :: b_incx, b_incy, b_n
! ......................................................................
!
    if ((nint(indm) .gt. 0.d0) .and. (nint(indp) .gt. 0.d0)) then
!
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigm, b_incx, tensm, b_incy)
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp, b_incx, tensp, b_incy)
!
        zernor = r8prem()
        if (icine .eq. 1) then
!
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, xm, b_incx, tensm, &
                       b_incy)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, xp, b_incx, tensp, &
                       b_incy)
!
        end if
!
!        PART HYDROSTATIQUE DES CONTRAINTES
        trt1 = (tensm(1)+tensm(2)+tensm(3))/3.d0
        trt2 = (tensp(1)+tensp(2)+tensp(3))/3.d0
!        PART DEVIATORIQUE DES CONTRAINTES
!
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tensm, b_incx, devm, b_incy)
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tensp, b_incx, devp, b_incy)
        do i = 1, 3
            devm(i) = devm(i)-trt1
            devp(i) = devp(i)-trt2
        end do
!         CALL DSCAL(NBSIG-3,SQRT(2.D0),DEVM(4),1)
!         CALL DSCAL(NBSIG-3,SQRT(2.D0),DEVP(4),1)
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        smeq = sqrt(ddot(b_n, devm, b_incx, devm, b_incy))
        b_n = to_blas_int(nbsig)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        speq = sqrt(ddot(b_n, devp, b_incx, devp, b_incy))
        preci = 1.d3*zernor*max(smeq, speq)
!
! ----      DANS LE CAS OU NORME(TENSM) = 0  OU NORME(DTENSMA) = 0 :
! ----      ON MET L'INDICATEUR A 0 :
!           -----------------------
        if (smeq .le. preci .or. speq .le. preci) then
            normdn = 0.d0
        else
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, devm, b_incx, n1, b_incy)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, devp, b_incx, n2, b_incy)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/smeq, n1, b_incx)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/speq, n2, b_incx)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, n1, b_incx, dn1n2, b_incy)
!
!          CALCUL DE LA DIFFERENCE N1 - N2
!
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, n2, b_incx, dn1n2, &
                       b_incy)
            b_n = to_blas_int(nbsig)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            ndn = sqrt(ddot(b_n, dn1n2, b_incx, dn1n2, b_incy))
            normdn = ndn/2.d0
        end if
!
    else
!        PAS D'INDICATEUR. ON NE FAIT RIEN
        normdn = 0.d0
!
    end if
!
end subroutine
