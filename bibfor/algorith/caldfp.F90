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
subroutine caldfp(msns, gamsns, dfpmdg, iret)
    implicit none
!
! person_in_charge: jean-michel.proix at edf.fr
!     ----------------------------------------------------------------
!
!     MONOCRISTAL : calcul des derivees de Fe en GDEF
!     IN  MSNS   : MS * NS
!         GAMSNS : somme des dGamma.Ms*Ns
!     OUT DFPMDG : dFp/dGamma_S
!         IRET   :  CODE RETOUR
!
!    Pour un seul systeme on calcule DF.FEn.d(Fp-1)/dGammaS
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/lcdetf.h"
#include "asterfort/matinv.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    integer(kind=8) :: iret, iopt, i, j, k, l
    real(kind=8) :: msns(3, 3), gamsns(3, 3), ddetdg
    real(kind=8) :: id(3, 3), coef, expo
    real(kind=8) :: det2, dfpmdg(3, 3)
    real(kind=8) :: dfpdg(3, 3), dfpmdf(3, 3, 3, 3), amax, amin, bmax, bmin
    real(kind=8) :: a(3, 3), am(3, 3), amt(3, 3), deta, coef2
    real(kind=8) :: b(3, 3), bm(3, 3), bmt(3, 3), detb
    blas_int :: b_incx, b_incy, b_n
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!     ----------------------------------------------------------------
!
    iret = 0
    iopt = 2
!
    if (iopt .eq. 1) then
!
!        calcul de dFp/dGamma suivant ANNAND 1996
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, gamsns, b_incx, a, b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, id, b_incx, a, &
                   b_incy)
!
!        TEST ANALOGUE A SIMO_MIEHE NMGPFI
        amax = 0.d0
        amin = 100.d0
        do i = 1, 3
            if (a(i, i) .gt. amax) amax = a(i, i)
            if (a(i, i) .lt. amin) amin = a(i, i)
        end do
        if ((amax .gt. 1.d3) .or. (amin .lt. 1.d-3)) then
            iret = 1
            goto 9999
        end if
!
        call lcdetf(3, a, deta)
!
        if (deta .gt. r8prem()) then
            expo = -1.d0/3.d0
            coef = deta**expo
        else
            iret = 1
            goto 9999
        end if
!
        call matinv('S', 3, a, am, det2)
        amt = transpose(am)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ddetdg = ddot(b_n, amt, b_incx, msns, b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, ddetdg, a, b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, a, b_incx, dfpdg, b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0/3.d0, dfpdg, b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, msns, b_incx, dfpdg, &
                   b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, coef, dfpdg, b_incx)
!
! calcul de dFp-1
        call r8inir(81, 0.d0, dfpmdf, 1)
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do l = 1, 3
                        dfpmdf(i, j, k, l) = dfpmdf(i, j, k, l)+am(i, k)*amt(j, l)
                    end do
                end do
            end do
        end do
        coef2 = -deta**(2.d0/3.d0)
!
        b_n = to_blas_int(81)
        b_incx = to_blas_int(1)
        call dscal(b_n, coef2, dfpmdf, b_incx)
!
        call r8inir(9, 0.d0, dfpmdg, 1)
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do l = 1, 3
                        dfpmdg(i, j) = dfpmdg(i, j)+dfpmdf(i, j, k, l)*dfpdg(k, l)
                    end do
                end do
            end do
        end do
!
    else if (iopt .eq. 2) then
!
!        calcul de dFp/dGamma par linearisation directe
!        de exp(-dgamma.ms x ns)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, gamsns, b_incx, b, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0, b, b_incx)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, id, b_incx, b, &
                   b_incy)
!
        bmax = 0.d0
        bmin = 100.d0
        do i = 1, 3
            if (b(i, i) .gt. bmax) bmax = b(i, i)
            if (b(i, i) .lt. bmin) bmin = b(i, i)
        end do
        if ((bmax .gt. 1.d3) .or. (bmin .lt. 1.d-3)) then
            iret = 1
            goto 9999
        end if
!
        call lcdetf(3, b, detb)
!
        if (detb .gt. r8prem()) then
            expo = -1.d0/3.d0
            coef = detb**expo
        else
            iret = 1
            goto 9999
        end if
!
        call matinv('S', 3, b, bm, det2)
!
        bmt = transpose(bm)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ddetdg = ddot(b_n, bmt, b_incx, msns, b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, ddetdg, b, b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, b, b_incx, dfpmdg, b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0/3.d0, dfpmdg, b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, msns, b_incx, dfpmdg, &
                   b_incy)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, -coef, dfpmdg, b_incx)
!
!
    else if (iopt .eq. 3) then
!
! suivant DE SOUZA-NIETO
        ASSERT(.false.)
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
!         CALL DEXPMAP(DFPMDG,NOCONV,DFP)
!
    else
        ASSERT(.false.)
!
    end if
!
!
9999 continue
end subroutine
