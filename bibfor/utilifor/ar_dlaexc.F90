! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
! ===============================================================
! THIS LAPACK 2.0 ROUTINE IS DEPRECATED
! DO NOT USE IT : YOU SHOULD PREFER UP-TO-DATE LAPACK ROUTINE
!
! BUT DO NOT REMOVE IT :
! THE PRESENT ROUTINE IS MANDATORY FOR ARPACK LIBRARY
! WHICH STICKS TO LAPACK 2.0 VERSION
! ==============================================================
subroutine ar_dlaexc(wantq, n, t, ldt, q, &
                     ldq, j1, n1, n2, work, &
                     info)
!
!     SUBROUTINE LAPACK PERMUTANT DEUX BLOCS DIAGONAUX D'UN MATRICE
!     TRIANGULAIRE SUPERIEURE.
!-----------------------------------------------------------------------
!  -- LAPACK AUXILIARY ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992
!
!  PURPOSE
!  =======
!
!  DLAEXC SWAPS ADJACENT DIAGONAL BLOCKS T11 AND T22 OF ORDER 1 OR 2 IN
!  AN UPPER QUASI-TRIANGULAR MATRIX T BY AN ORTHOGONAL SIMILARITY
!  TRANSFORMATION.
!
!  T MUST BE IN SCHUR CANONICAL FORM, THAT IS, BLOCK UPPER TRIANGULAR
!  WITH 1-BY-1 AND 2-BY-2 DIAGONAL BLOCKS, EACH 2-BY-2 DIAGONAL BLOCK
!  HAS ITS DIAGONAL ELEMNTS EQUAL AND ITS OFF-DIAGONAL ELEMENTS OF
!  OPPOSITE SIGN.
!
!  ARGUMENTS
!  =========
!
!  WANTQ   (INPUT) LOGICAL
!          = .TRUE. : ACCUMULATE THE TRANSFORMATION IN THE MATRIX Q,
!          = .FALSE.: DO NOT ACCUMULATE THE TRANSFORMATION.
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX T. N >= 0.
!
!  T       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDT,N)
!          ON ENTRY, THE UPPER QUASI-TRIANGULAR MATRIX T, IN SCHUR
!          CANONICAL FORM.
!          ON EXIT, THE UPDATED MATRIX T, AGAIN IN SCHUR CANONICAL FORM.
!
!  LDT     (INPUT)  INTEGER
!          THE LEADING DIMENSION OF THE ARRAY T. LDT >= MAX(1,N).
!
!  Q       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDQ,N)
!          ON ENTRY, IF WANTQ IS .TRUE., THE ORTHOGONAL MATRIX Q.
!          ON EXIT, IF WANTQ IS .TRUE., THE UPDATED MATRIX Q.
!          IF WANTQ IS .FALSE., Q IS NOT REFERENCED.
!
!  LDQ     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY Q.
!          LDQ >= 1, AND IF WANTQ IS .TRUE., LDQ >= N.
!
!  J1      (INPUT) INTEGER
!          THE INDEX OF THE FIRST ROW OF THE FIRST BLOCK T11.
!
!  N1      (INPUT) INTEGER
!          THE ORDER OF THE FIRST BLOCK T11. N1 = 0, 1 OR 2.
!
!  N2      (INPUT) INTEGER
!          THE ORDER OF THE SECOND BLOCK T22. N2 = 0, 1 OR 2.
!
!  WORK    (WORKSPACE) REAL*8 ARRAY, DIMENSION (N)
!
!  INFO    (OUTPUT) INTEGER
!          = 0: SUCCESSFUL EXIT
!          = 1: THE TRANSFORMED MATRIX T WOULD BE TOO FAR FROM SCHUR
!               FORM, THE BLOCKS ARE NOT SWAPPED AND T AND Q ARE
!               UNCHANGED.
!
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 3 RETURN PAR GOTO 999,
!            REMPLACEMENT DE DLAMCH PAR R8PREM, R8MIEM ET ISBAEM,
!            MODIFICATION DES APPELS BLAS (ROUTINE ASTER BL...),
!            IMPLICIT NONE.
! INTRINSIC FUNCTIONS
!            ABS, MAX.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/ar_dlanv2.h"
#include "asterfort/ar_dlarfg.h"
#include "asterfort/ar_dlartg.h"
#include "asterfort/ar_dlasy2.h"
#include "blas/dlacpy.h"
#include "blas/dlange.h"
#include "blas/dlarfx.h"
#include "blas/drot.h"
    aster_logical :: wantq
    integer(kind=8) :: info, j1, ldq, ldt, n, n1, n2
!     ..
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: q(ldq, *), t(ldt, *), work(*)
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: zero, one
    parameter(zero=0.0d+0, one=1.0d+0)
    real(kind=8) :: ten
    parameter(ten=1.0d+1)
    integer(kind=8) :: ldd, ldx
    parameter(ldd=4, ldx=2)
!     ..
!     .. LOCAL SCALARS ..
    integer(kind=8) :: ierr, j2, j3, j4, k, nd
    real(kind=8) :: cs, dnorm, eps, scale, smlnum, sn, t11, t22, t33, tau, tau1
    real(kind=8) :: tau2, temp, thresh, wi1, wi2, wr1, wr2, xnorm
!     ..
!     .. LOCAL ARRAYS ..
    real(kind=8) :: d(ldd, 4), u(3), u1(3), u2(3), x(ldx, 2)
    blas_int :: b_lda, b_ldb, b_m, b_n
    blas_int :: b_ldc
    blas_int :: b_incx, b_incy
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
    info = 0
!
!     QUICK RETURN IF POSSIBLE
!
    if (n .eq. 0 .or. n1 .eq. 0 .or. n2 .eq. 0) goto 999
    if (j1+n1 .gt. n) goto 999
!
    j2 = j1+1
    j3 = j1+2
    j4 = j1+3
!
    if (n1 .eq. 1 .and. n2 .eq. 1) then
!
!        SWAP TWO 1-BY-1 BLOCKS.
!
        t11 = t(j1, j1)
        t22 = t(j2, j2)
!
!        DETERMINE THE TRANSFORMATION TO PERFORM THE INTERCHANGE.
!
        call ar_dlartg(t(j1, j2), t22-t11, cs, sn, temp)
!
!        APPLY TRANSFORMATION TO THE MATRIX T.
!
        if (j3 .le. n) then
            b_n = to_blas_int(n-j1-1)
            b_incx = to_blas_int(ldt)
            b_incy = to_blas_int(ldt)
            call drot(b_n, t(j1, j3), b_incx, t(j2, j3), b_incy, &
                      cs, sn)
        end if
        b_n = to_blas_int(j1-1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call drot(b_n, t(1, j1), b_incx, t(1, j2), b_incy, &
                  cs, sn)
!
        t(j1, j1) = t22
        t(j2, j2) = t11
!
        if (wantq) then
!
!           ACCUMULATE TRANSFORMATION IN THE MATRIX Q.
!
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call drot(b_n, q(1, j1), b_incx, q(1, j2), b_incy, &
                      cs, sn)
        end if
!
    else
!
!        SWAPPING INVOLVES AT LEAST ONE 2-BY-2 BLOCK.
!
!        COPY THE DIAGONAL BLOCK OF ORDER N1+N2 TO THE LOCAL ARRAY D
!        AND COMPUTE ITS NORM.
!
        nd = n1+n2
! DUE TO CRP102 CALL DLACPY( 'FULL', ND, ND, T( J1, J1 ), LDT, D, LDD )
        b_ldb = to_blas_int(ldd)
        b_lda = to_blas_int(ldt)
        b_m = to_blas_int(nd)
        b_n = to_blas_int(nd)
        call dlacpy('F', b_m, b_n, t(j1, j1), b_lda, &
                    d, b_ldb)
! DUE TO CRP102 DNORM = DLANGE( 'MAX', ND, ND, D, LDD, WORK )
        b_lda = to_blas_int(ldd)
        b_m = to_blas_int(nd)
        b_n = to_blas_int(nd)
        dnorm = dlange('M', b_m, b_n, d, b_lda, work)
!
!        COMPUTE MACHINE-DEPENDENT THRESHOLD FOR TEST FOR ACCEPTING
!        SWAP.
!
        eps = r8prem()*0.5d0*isbaem()
        smlnum = r8miem()/eps
        thresh = max(ten*eps*dnorm, smlnum)
!
!        SOLVE T11*X - X*T22 = SCALE*T12 FOR X.
!
        call ar_dlasy2(.false._1, .false._1, -1, n1, n2, &
                       d, ldd, d(n1+1, n1+1), ldd, d(1, n1+1), &
                       ldd, scale, x, ldx, xnorm, &
                       ierr)
!
!        SWAP THE ADJACENT DIAGONAL BLOCKS.
!
        k = n1+n1+n2-3
!
        select case (k)
        case (1)
!
!           N1 = 1, N2 = 2: GENERATE ELEMENTARY REFLECTOR H SO THAT:
!
!           ( SCALE, X11, X12 ) H = ( 0, 0, * )
!
            u(1) = scale
            u(2) = x(1, 1)
            u(3) = x(1, 2)
            call ar_dlarfg(3, u(3), u, 1, tau)
            u(3) = one
            t11 = t(j1, j1)
!
!           PERFORM SWAP PROVISIONALLY ON DIAGONAL BLOCK IN D.
!
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            call dlarfx('L', b_m, b_n, u, tau, &
                        d, b_ldc, work)
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u, tau, &
                        d, b_ldc, work)
!
!           TEST WHETHER TO REJECT SWAP.
!
            if (max(abs(d(3, 1)), abs(d(3, 2)), abs(d(3, 3)-t11)) .gt. thresh) goto 50
!
!           ACCEPT SWAP: APPLY TRANSFORMATION TO THE ENTIRE MATRIX T.
!
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(3)
            b_n = to_blas_int(n-j1+1)
            call dlarfx('L', b_m, b_n, u, tau, &
                        t(j1, j1), b_ldc, work)
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(j2)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u, tau, &
                        t(1, j1), b_ldc, work)
!
            t(j3, j1) = zero
            t(j3, j2) = zero
            t(j3, j3) = t11
!
            if (wantq) then
!
!               ACCUMULATE TRANSFORMATION IN THE MATRIX Q.
!
                b_ldc = to_blas_int(ldq)
                b_m = to_blas_int(n)
                b_n = to_blas_int(3)
                call dlarfx('R', b_m, b_n, u, tau, &
                            q(1, j1), b_ldc, work)
            end if
!
        case (2)
!           N1 = 2, N2 = 1: GENERATE ELEMENTARY REFLECTOR H SO THAT:
!
!           H (  -X11 ) = ( * )
!             (  -X21 ) = ( 0 )
!             ( SCALE ) = ( 0 )
!
            u(1) = -x(1, 1)
            u(2) = -x(2, 1)
            u(3) = scale
            call ar_dlarfg(3, u(1), u(2), 1, tau)
            u(1) = one
            t33 = t(j3, j3)
!
!            PERFORM SWAP PROVISIONALLY ON DIAGONAL BLOCK IN D.
!
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            call dlarfx('L', b_m, b_n, u, tau, &
                        d, b_ldc, work)
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u, tau, &
                        d, b_ldc, work)
!
!            TEST WHETHER TO REJECT SWAP.
!
            if (max(abs(d(2, 1)), abs(d(3, 1)), abs(d(1, 1)-t33)) .gt. thresh) goto 50
!
!            ACCEPT SWAP: APPLY TRANSFORMATION TO THE ENTIRE MATRIX T.
!
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(j3)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u, tau, &
                        t(1, j1), b_ldc, work)
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(3)
            b_n = to_blas_int(n-j1)
            call dlarfx('L', b_m, b_n, u, tau, &
                        t(j1, j2), b_ldc, work)
!
            t(j1, j1) = t33
            t(j2, j1) = zero
            t(j3, j1) = zero
!
            if (wantq) then
!
!               ACCUMULATE TRANSFORMATION IN THE MATRIX Q.
!
                b_ldc = to_blas_int(ldq)
                b_m = to_blas_int(n)
                b_n = to_blas_int(3)
                call dlarfx('R', b_m, b_n, u, tau, &
                            q(1, j1), b_ldc, work)
            end if
!
        case (3)
!
!            N1 = 2, N2 = 2: GENERATE ELEMENTARY REFLECTORS H(1) AND H(2) SO
!            THAT:
!
!            H(2) H(1) (  -X11  -X12 ) = (  *  * )
!                      (  -X21  -X22 )   (  0  * )
!                      ( SCALE    0  )   (  0  0 )
!                      (    0  SCALE )   (  0  0 )
!
            u1(1) = -x(1, 1)
            u1(2) = -x(2, 1)
            u1(3) = scale
            call ar_dlarfg(3, u1(1), u1(2), 1, tau1)
            u1(1) = one
!
            temp = -tau1*(x(1, 2)+u1(2)*x(2, 2))
            u2(1) = -temp*u1(2)-x(2, 2)
            u2(2) = -temp*u1(3)
            u2(3) = scale
            call ar_dlarfg(3, u2(1), u2(2), 1, tau2)
            u2(1) = one
!
!            PERFORM SWAP PROVISIONALLY ON DIAGONAL BLOCK IN D.
!
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(4)
            call dlarfx('L', b_m, b_n, u1, tau1, &
                        d, b_ldc, work)
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(4)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u1, tau1, &
                        d, b_ldc, work)
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(3)
            b_n = to_blas_int(4)
            call dlarfx('L', b_m, b_n, u2, tau2, &
                        d(2, 1), b_ldc, work)
            b_ldc = to_blas_int(ldd)
            b_m = to_blas_int(4)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u2, tau2, &
                        d(1, 2), b_ldc, work)
!
!            TEST WHETHER TO REJECT SWAP.
!
            if (max(abs(d(3, 1)), abs(d(3, 2)), abs(d(4, 1)), abs(d(4, 2))) .gt. thresh) &
                goto 50
!
!            ACCEPT SWAP: APPLY TRANSFORMATION TO THE ENTIRE MATRIX T.
!
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(3)
            b_n = to_blas_int(n-j1+1)
            call dlarfx('L', b_m, b_n, u1, tau1, &
                        t(j1, j1), b_ldc, work)
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(j4)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u1, tau1, &
                        t(1, j1), b_ldc, work)
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(3)
            b_n = to_blas_int(n-j1+1)
            call dlarfx('L', b_m, b_n, u2, tau2, &
                        t(j2, j1), b_ldc, work)
            b_ldc = to_blas_int(ldt)
            b_m = to_blas_int(j4)
            b_n = to_blas_int(3)
            call dlarfx('R', b_m, b_n, u2, tau2, &
                        t(1, j2), b_ldc, work)
!
            t(j3, j1) = zero
            t(j3, j2) = zero
            t(j4, j1) = zero
            t(j4, j2) = zero
!
            if (wantq) then
!
!               ACCUMULATE TRANSFORMATION IN THE MATRIX Q.
!
                b_ldc = to_blas_int(ldq)
                b_m = to_blas_int(n)
                b_n = to_blas_int(3)
                call dlarfx('R', b_m, b_n, u1, tau1, &
                            q(1, j1), b_ldc, work)
                b_ldc = to_blas_int(ldq)
                b_m = to_blas_int(n)
                b_n = to_blas_int(3)
                call dlarfx('R', b_m, b_n, u2, tau2, &
                            q(1, j2), b_ldc, work)
            end if
!
        end select
!
        if (n2 .eq. 2) then
!
!           STANDARDIZE NEW 2-BY-2 BLOCK T11
!
            call ar_dlanv2(t(j1, j1), t(j1, j2), t(j2, j1), t(j2, j2), wr1, &
                           wi1, wr2, wi2, cs, sn)
            b_n = to_blas_int(n-j1-1)
            b_incx = to_blas_int(ldt)
            b_incy = to_blas_int(ldt)
            call drot(b_n, t(j1, j1+2), b_incx, t(j2, j1+2), b_incy, &
                      cs, sn)
            b_n = to_blas_int(j1-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call drot(b_n, t(1, j1), b_incx, t(1, j2), b_incy, &
                      cs, sn)
            if (wantq) then
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call drot(b_n, q(1, j1), b_incx, q(1, j2), b_incy, &
                          cs, sn)
            end if
        end if
!
        if (n1 .eq. 2) then
!
!           STANDARDIZE NEW 2-BY-2 BLOCK T22
!
            j3 = j1+n2
            j4 = j3+1
            call ar_dlanv2(t(j3, j3), t(j3, j4), t(j4, j3), t(j4, j4), wr1, &
                           wi1, wr2, wi2, cs, sn)
            if (j3+2 .le. n) then
                b_n = to_blas_int(n-j3-1)
                b_incx = to_blas_int(ldt)
                b_incy = to_blas_int(ldt)
                call drot(b_n, t(j3, j3+2), b_incx, t(j4, j3+2), b_incy, &
                          cs, sn)
            end if
            b_n = to_blas_int(j3-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call drot(b_n, t(1, j3), b_incx, t(1, j4), b_incy, &
                      cs, sn)
            if (wantq) then
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call drot(b_n, q(1, j3), b_incx, q(1, j4), b_incy, &
                          cs, sn)
            end if
        end if
!
    end if
    goto 999
!
!     EXIT WITH INFO = 1 IF SWAP WAS REJECTED.
!
50  continue
    info = 1
999 continue
!
    call matfpe(1)
!
!     END OF DLAEXC
!
end subroutine
