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
subroutine calcfe(nr, ndt, nvi, vind, df, &
                  gamsns, fe, fp, iret)
    implicit none
!       MONOCRISTAL : CALCUL DE Fe et Fp, F=Fe.Fp
!       IN  NR     :  DIMENSION DECLAREE DRDY
!           NDT    :  NOMBRE DE COMPOSANTES DE SIGMA (6)
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!                     contiennent Fe(t), Fp(t)
!           DF     :  Increment de Gradient de deformation
!           GAMSNS :  Somme de dGamma.Ms*Ns
!       OUT FE     :  Gradient de tranformation elastique
!           FP     :  Gradient de tranformation plastique
!
#include "asterc/r8prem.h"
#include "asterfort/lcdetf.h"
#include "asterfort/matinv.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    real(kind=8) :: fe(3, 3), df(3, 3), gamsns(3, 3), dffe(3, 3), dfpm(3, 3)
    real(kind=8) :: fem(3, 3)
    real(kind=8) :: vind(*), id(3, 3), det, coef, dfp(3, 3), expo, fp(3, 3)
    real(kind=8) :: fpm(3, 3), dfpmax, dfpmin, det2
    integer(kind=8) :: nr, ndt, iret, iopt, i, nvi
    blas_int :: b_incx, b_incy, b_n
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!     ----------------------------------------------------------------
!
    iret = 0
    iopt = 2
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vind(nvi-3-18+10), b_incx, fem, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vind(nvi-3-18+1), b_incx, fpm, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, id, b_incx, fem, &
               b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, id, b_incx, fpm, &
               b_incy)
!
    dffe = matmul(df, fem)
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, gamsns, b_incx, dfp, b_incy)
!
    if (iopt .eq. 1) then
!
!        suivant ANNAND 1996
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, id, b_incx, dfp, &
                   b_incy)
!
!        TEST ANALOGUE A SIMO_MIEHE NMGPFI
        dfpmax = 0.d0
        dfpmin = 100.d0
        do i = 1, 3
            if (dfp(i, i) .gt. dfpmax) dfpmax = dfp(i, i)
            if (dfp(i, i) .lt. dfpmin) dfpmin = dfp(i, i)
        end do
        if ((dfpmax .gt. 1.d3) .or. (dfpmin .lt. 1.d-3)) then
            iret = 1
            goto 999
        end if
!
        call lcdetf(3, dfp, det)
!
        if (det .gt. r8prem()) then
            expo = -1.d0/3.d0
            coef = det**expo
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            call dscal(b_n, coef, dfp, b_incx)
        else
            iret = 1
            goto 999
        end if
!
        call matinv('S', 3, dfp, dfpm, det2)
!
    else if (iopt .eq. 2) then
!
!        linearisation directe de exp(-dgamma.ms x ns)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dfp, b_incx, dfpm, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0, dfpm, b_incx)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, id, b_incx, dfpm, &
                   b_incy)
!
        dfpmax = 0.d0
        dfpmin = 100.d0
        do i = 1, 3
            if (dfpm(i, i) .gt. dfpmax) dfpmax = dfpm(i, i)
            if (dfpm(i, i) .lt. dfpmin) dfpmin = dfpm(i, i)
        end do
        if ((dfpmax .gt. 1.d3) .or. (dfpmin .lt. 1.d-3)) then
            iret = 1
            goto 999
        end if
!
        call lcdetf(3, dfpm, det)
!
        if (det .gt. r8prem()) then
            expo = -1.d0/3.d0
            coef = det**expo
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            call dscal(b_n, coef, dfpm, b_incx)
        else
            iret = 1
            goto 999
        end if
!
        call matinv('S', 3, dfpm, dfp, det2)
!
!
    else if (iopt .eq. 3) then
!
! suivant DE SOUZA-NIETO
!
!         DFPMAX=0.D0
!         DFPMIN=100.D0
!         DO 30 I=1,3
!            IF (DFP(I,I).GT.DFPMAX) DFPMAX=DFP(I,I)
!            IF (DFP(I,I).LT.DFPMIN) DFPMIN=DFP(I,I)
! 30      CONTINUE
!         IF ((ABS(DFPMAX).GT.10.D0).OR.(ABS(DFPMIN).GT.10.D0)) THEN
!           IRET=1
!           GOTO 9999
!         ENDIF
!         CALL DSCAL(9,-1.0D0,DFP,1)
!         CALL EXPMAP(DFPM,NOCONV,DFP)
!         IF (NOCONV) THEN
!            IRET=1
!            GOTO 9999
!         ENDIF
!         CALL MATINV('S',3,DFPM,DFP,DET2)
!
    end if
!
    fe = matmul(dffe, dfpm)
!
! post traitement
!
    fp = matmul(dfp, fpm)
!
999 continue
!
end subroutine
