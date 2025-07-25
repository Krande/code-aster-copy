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
subroutine znapps(n, kev, np, shift, v, &
                  ldv, h, ldh, resid, q, &
                  ldq, workl, workd)
    implicit none
!
!     SUBROUTINE ARPACK PREPARANT LE RESTART VIA UN QR IMPLICITE POUR
!     ELIMINER LES NP MODES PROPRES INDESIRABLES.
!-----------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZNAPPS
!
!\DESCRIPTION:
!  GIVEN THE ARNOLDI FACTORIZATION
!
!     A*V_{K} - V_{K}*H_{K} = R_{K+P}*E_{K+P}^T,
!
!  APPLY NP IMPLICIT SHIFTS RESULTING IN
!
!     A*(V_{K}*Q) - (V_{K}*Q)*(Q^T* H_{K}*Q) = R_{K+P}*E_{K+P}^T * Q
!
!  WHERE Q IS AN ORTHOGONAL MATRIX WHICH IS THE PRODUCT OF ROTATIONS
!  AND REFLECTIONS RESULTING FROM THE NP BULGE CHANGE SWEEPS.
!  THE UPDATED ARNOLDI FACTORIZATION BECOMES:
!
!     A*VNEW_{K} - VNEW_{K}*HNEW_{K} = RNEW_{K}*E_{K}^T.
!
!\USAGE:
!  CALL ZNAPPS
!     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ,
!       WORKL, WORKD )
!
!\ARGUMENTS
!  N       INTEGER.  (INPUT)
!          PROBLEM SIZE, I.E. SIZE OF MATRIX A.
!
!  KEV     INTEGER.  (INPUT/OUTPUT)
!          KEV+NP IS THE SIZE OF THE INPUT MATRIX H.
!          KEV IS THE SIZE OF THE UPDATED MATRIX HNEW.
!
!  NP      INTEGER.  (INPUT)
!          NUMBER OF IMPLICIT SHIFTS TO BE APPLIED.
!
!  SHIFT   COMPLEX*16 ARRAY OF LENGTH NP.  (INPUT)
!          THE SHIFTS TO BE APPLIED.
!
!  V       COMPLEX*16 N BY (KEV+NP) ARRAY.  (INPUT/OUTPUT)
!          ON INPUT, V CONTAINS THE CURRENT KEV+NP ARNOLDI VECTORS.
!          ON OUTPUT, V CONTAINS THE UPDATED KEV ARNOLDI VECTORS
!          IN THE FIRST KEV COLUMNS OF V.
!
!  LDV     INTEGER.  (INPUT)
!          LEADING DIMENSION OF V EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  H       COMPLEX*16 (KEV+NP) BY (KEV+NP) ARRAY.  (INPUT/OUTPUT)
!          ON INPUT, H CONTAINS THE CURRENT KEV+NP BY KEV+NP UPPER
!          HESSENBERG MATRIX OF THE ARNOLDI FACTORIZATION.
!          ON OUTPUT, H CONTAINS THE UPDATED KEV BY KEV UPPER HESSENBERG
!          MATRIX IN THE KEV LEADING SUBMATRIX.
!
!  LDH     INTEGER.  (INPUT)
!          LEADING DIMENSION OF H EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  RESID   COMPLEX*16 ARRAY OF LENGTH N.  (INPUT/OUTPUT)
!          ON INPUT, RESID CONTAINS THE THE RESIDUAL VECTOR R_{K+P}.
!          ON OUTPUT, RESID IS THE UPDATE RESIDUAL VECTOR RNEW_{K}
!          IN THE FIRST KEV LOCATIONS.
!
!  Q       COMPLEX*16 KEV+NP BY KEV+NP WORK ARRAY.  (WORKSPACE)
!          WORK ARRAY USED TO ACCUMULATE THE ROTATIONS AND REFLECTIONS
!          DURING THE BULGE CHASE SWEEP.
!
!  LDQ     INTEGER.  (INPUT)
!          LEADING DIMENSION OF Q EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  WORKL   COMPLEX*16 WORK ARRAY OF LENGTH (KEV+NP).  (WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.
!
!  WORKD   COMPLEX*16 WORK ARRAY OF LENGTH 2*N.  (WORKSPACE)
!          DISTRIBUTED ARRAY USED IN THE APPLICATION OF THE ACCUMULATED
!          ORTHOGONAL MATRIX Q.
!
!\ENDDOC
!
!-----------------------------------------------------------------------
!
!\BEGINLIB
!
!\LOCAL VARIABLES:
!     XXXXXX  COMPLEX*16
!
!\REFERENCES:
!  1. D.C. SORENSEN, "IMPLICIT APPLICATION OF POLYNOMIAL FILTERS IN
!     A K-STEP ARNOLDI METHOD", SIAM J. MATR. ANAL. APPS., 13 (1992),
!     PP 357-385.
!
!\ROUTINES CALLED:
!     IVOUT   ARPACK UTILITY ROUTINE THAT PRINTS INTEGERS.
!     SECOND  ARPACK UTILITY ROUTINE FOR TIMING.
!     ZMOUT   ARPACK UTILITY ROUTINE THAT PRINTS MATRICES
!     ZVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     ZLACPY  LAPACK MATRIX COPY ROUTINE.
!     ZLANHS  LAPACK ROUTINE THAT COMPUTES VARIOUS NORMS OF A MATRIX.
!     ZLARTG  LAPACK GIVENS ROTATION CONSTRUCTION ROUTINE.
!     ZLASET  LAPACK MATRIX INITIALIZATION ROUTINE.
!     DLABAD  LAPACK ROUTINE FOR DEFINING THE UNDERFLOW AND OVERFLOW
!             LIMITS.
!     DLAMCH  LAPACK ROUTINE THAT DETERMINES MACHINE CONSTANTS.
!     DLAPY2  LAPACK ROUTINE TO COMPUTE SQRT(X**2+Y**2) CAREFULLY.
!     ZGEMV   LEVEL 2 BLAS ROUTINE FOR MATRIX VECTOR MULTIPLICATION.
!     ZAXPY   LEVEL 1 BLAS THAT COMPUTES A VECTOR TRIAD.
!     ZCOPY   LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER.
!     ZLSCAL   LEVEL 1 BLAS THAT SCALES A VECTOR.
!
!\AUTHOR
!     DANNY SORENSEN               PHUONG VU
!     RICHARD LEHOUCQ              CRPC / RICE UNIVERSITY
!     DEPT. OF COMPUTATIONAL &     HOUSTON, TEXAS
!     APPLIED MATHEMATICS
!     RICE UNIVERSITY
!     HOUSTON, TEXAS
!
!\SCCS INFORMATION: @(#)
! FILE: NAPPS.F   SID: 2.3   DATE OF SID: 3/28/97   RELEASE: 2
!
!\REMARKS
!  1. IN THIS VERSION, EACH SHIFT IS APPLIED TO ALL THE SUBLOCKS OF
!     THE HESSENBERG MATRIX H AND NOT JUST TO THE SUBMATRIX THAT IT
!     COMES FROM. DEFLATION AS IN LAPACK ROUTINE ZLAHQR (QR ALGORITHM
!     FOR UPPER HESSENBERG MATRICES ) IS USED.
!     UPON OUTPUT, THE SUBDIAGONALS OF H ARE ENFORCED TO BE NON-NEGATIVE
!     REAL NUMBERS.
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
!
!     %-----------------------------%
!     | INCLUDE FILES FOR DEBUGGING |
!     %-----------------------------%
!
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/ar_zlartg.h"
#include "asterfort/ivout.h"
#include "asterfort/zlscal.h"
#include "asterfort/zmout.h"
#include "asterfort/zvout.h"
#include "blas/dlapy2.h"
#include "blas/zaxpy.h"
#include "blas/zcopy.h"
#include "blas/zgemv.h"
#include "blas/zlacpy.h"
#include "blas/zlanhs.h"
#include "blas/zlaset.h"
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr, mneigh, mnapps
    integer(kind=8) :: mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
    integer(kind=8) :: kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    complex(kind=8) :: h(ldh, kev+np), resid(n), shift(np), v(ldv, kev+np)
    complex(kind=8) :: q(ldq, kev+np), workd(2*n), workl(kev+np)
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    complex(kind=8) :: one, zero
    real(kind=8) :: rzero
    parameter(one=(1.0d+0, 0.0d+0), zero=(0.0d+0, 0.0d+0),&
     &           rzero=0.0d+0)
!
!     %------------------------%
!     | LOCAL SCALARS & ARRAYS |
!     %------------------------%
!
    integer(kind=8) :: i, iend, istart, j, jj, kplusp, msglvl
    aster_logical :: first
    complex(kind=8) :: f, g, h11, h21, r, s, sigma(1), t
    real(kind=8) :: c, smlnum, ulp, unfl, tst1, rbid(1)
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda, b_m
    blas_int :: b_ldb
    save first, smlnum, ulp, unfl
!
!
!     %--------------------%
!     | EXTERNAL FUNCTIONS |
!     %--------------------%
!
!
!     %---------------------%
!     | STATEMENT FUNCTIONS |
!     %---------------------%
!
#define zabs1( cdum ) abs( dble( cdum ) ) + abs( dimag( cdum ) )
!
!     %----------------%
!     | DATA STATMENTS |
!     %----------------%
!
    data first/.true./
!
!     %-----------------------%
!     | EXECUTABLE STATEMENTS |
!     %-----------------------%
!
    call matfpe(-1)
!
    rbid = 0.d0
    if (first) then
!
!        %-----------------------------------------------%
!        | SET MACHINE-DEPENDENT CONSTANTS FOR THE       |
!        | STOPPING CRITERION. IF NORM(H) <= SQRT(OVFL), |
!        | OVERFLOW SHOULD NOT OCCUR.                    |
!        | REFERENCE: LAPACK SUBROUTINE ZLAHQR           |
!        %-----------------------------------------------%
!
!
        unfl = r8miem()
! DUE RO CRS512         OVFL = ONE / UNFL
        ulp = r8prem()*0.5d0*isbaem()
        smlnum = unfl*(n/ulp)
        first = .false.
    end if
!
!     %-------------------------------%
!     | INITIALIZE TIMING STATISTICS  |
!     | & MESSAGE LEVEL FOR DEBUGGING |
!     %-------------------------------%
!
    msglvl = mnapps
!
    kplusp = kev+np
!
!     %--------------------------------------------%
!     | INITIALIZE Q TO THE IDENTITY TO ACCUMULATE |
!     | THE ROTATIONS AND REFLECTIONS              |
!     %--------------------------------------------%
!
    b_lda = to_blas_int(ldq)
    b_m = to_blas_int(kplusp)
    b_n = to_blas_int(kplusp)
    call zlaset('A', b_m, b_n, zero, one, &
                q, b_lda)
!
!     %----------------------------------------------%
!     | QUICK RETURN IF THERE ARE NO SHIFTS TO APPLY |
!     %----------------------------------------------%
!
    if (np .eq. 0) goto 9000
!
!     %----------------------------------------------%
!     | CHASE THE BULGE WITH THE APPLICATION OF EACH |
!     | IMPLICIT SHIFT. EACH SHIFT IS APPLIED TO THE |
!     | WHOLE MATRIX INCLUDING EACH BLOCK.           |
!     %----------------------------------------------%
!
    do jj = 1, np
        sigma(1) = shift(jj)
!
        if (msglvl .gt. 2) then
            call ivout(logfil, 1, [jj], ndigit, '_NAPPS: SHIFT NUMBER.')
            call zvout(logfil, 1, sigma, ndigit, '_NAPPS: VALUE OF THE SHIFT ')
        end if
!
        istart = 1
20      continue
!
        do i = istart, kplusp-1
!
!           %----------------------------------------%
!           | CHECK FOR SPLITTING AND DEFLATION. USE |
!           | A STANDARD TEST AS IN THE QR ALGORITHM |
!           | REFERENCE: LAPACK SUBROUTINE ZLAHQR    |
!           %----------------------------------------%
!
            tst1 = zabs1(h(i, i))+zabs1(h(i+1, i+1))
            b_lda = to_blas_int(ldh)
            b_n = to_blas_int(kplusp-jj+1)
            if (tst1 .eq. rzero) tst1 = zlanhs('1', b_n, h, b_lda, rbid)
!
            if (abs(dble(h(i+1, i))) .le. max(ulp*tst1, smlnum)) then
                if (msglvl .gt. 0) then
                    call ivout(logfil, 1, [i], ndigit, &
                               '_NAPPS: MATRIX SPLITTING AT ROW/COLUMN NO.')
                    call ivout(logfil, 1, [jj], ndigit, &
                               '_NAPPS: MATRIX SPLITTING WITH SHIFT NUMBER.')
                    call zvout(logfil, 1, h(i+1, i), ndigit, '_NAPPS: OFF DIAGONAL ELEMENT.')
                end if
                iend = i
                h(i+1, i) = zero
                goto 40
            end if
        end do
        iend = kplusp
40      continue
!
        if (msglvl .gt. 2) then
            call ivout(logfil, 1, [istart], ndigit, '_NAPPS: START OF CURRENT BLOCK ')
            call ivout(logfil, 1, [iend], ndigit, '_NAPPS: END OF CURRENT BLOCK ')
        end if
!
!        %------------------------------------------------%
!        | NO REASON TO APPLY A SHIFT TO BLOCK OF ORDER 1 |
!        | OR IF THE CURRENT BLOCK STARTS AFTER THE POINT |
!        | OF COMPRESSION SINCE WE'LL DISCARD THIS STUFF  |
!        %------------------------------------------------%
!
        if (istart .eq. iend .or. istart .gt. kev) goto 100
!
        h11 = h(istart, istart)
        h21 = h(istart+1, istart)
        f = h11-sigma(1)
        g = h21
!
        do i = istart, iend-1
!
!           %------------------------------------------------------%
!           | CONSTRUCT THE PLANE ROTATION G TO ZERO OUT THE BULGE |
!           %------------------------------------------------------%
!
            call ar_zlartg(f, g, c, s, r)
            if (i .gt. istart) then
                h(i, i-1) = r
                h(i+1, i-1) = zero
            end if
!
!           %---------------------------------------------%
!           | APPLY ROTATION TO THE LEFT OF H;  H <- G'*H |
!           %---------------------------------------------%
!
            do j = i, kplusp
                t = c*h(i, j)+s*h(i+1, j)
                h(i+1, j) = -dconjg(s)*h(i, j)+c*h(i+1, j)
                h(i, j) = t
            end do
!
!           %---------------------------------------------%
!           | APPLY ROTATION TO THE RIGHT OF H;  H <- H*G |
!           %---------------------------------------------%
!
            do j = 1, min(i+2, iend)
                t = c*h(j, i)+dconjg(s)*h(j, i+1)
                h(j, i+1) = -s*h(j, i)+c*h(j, i+1)
                h(j, i) = t
            end do
!
!           %-----------------------------------------------------%
!           | ACCUMULATE THE ROTATION IN THE MATRIX Q;  Q <- Q*G' |
!           %-----------------------------------------------------%
!
            do j = 1, min(i+jj, kplusp)
                t = c*q(j, i)+dconjg(s)*q(j, i+1)
                q(j, i+1) = -s*q(j, i)+c*q(j, i+1)
                q(j, i) = t
            end do
!
!           %---------------------------%
!           | PREPARE FOR NEXT ROTATION |
!           %---------------------------%
!
            if (i .lt. iend-1) then
                f = h(i+1, i)
                g = h(i+2, i)
            end if
        end do
!
!        %-------------------------------%
!        | FINISHED APPLYING THE SHIFT.  |
!        %-------------------------------%
!
100     continue
!
!        %---------------------------------------------------------%
!        | APPLY THE SAME SHIFT TO THE NEXT BLOCK IF THERE IS ANY. |
!        %---------------------------------------------------------%
!
        istart = iend+1
        if (iend .lt. kplusp) goto 20
!
!        %---------------------------------------------%
!        | LOOP BACK TO THE TOP TO GET THE NEXT SHIFT. |
!        %---------------------------------------------%
!
    end do
!
!     %---------------------------------------------------%
!     | PERFORM A SIMILARITY TRANSFORMATION THAT MAKES    |
!     | SURE THAT THE COMPRESSED H WILL HAVE NON-NEGATIVE |
!     | REAL SUBDIAGONAL ELEMENTS.                        |
!     %---------------------------------------------------%
!
    do j = 1, kev
        if (dble(h(j+1, j)) .lt. rzero .or. dimag(h(j+1, j)) .ne. rzero) then
            t = h(j+1, j)/dlapy2(dble(h(j+1, j)), dimag(h(j+1, j)))
            call zlscal(kplusp-j+1, dconjg(t), h(j+1, j), ldh)
            call zlscal(min(j+2, kplusp), t, h(1, j+1), 1)
            call zlscal(min(j+np+1, kplusp), t, q(1, j+1), 1)
            h(j+1, j) = dcmplx(dble(h(j+1, j)), rzero)
        end if
    end do
!
    do i = 1, kev
!
!        %--------------------------------------------%
!        | FINAL CHECK FOR SPLITTING AND DEFLATION.   |
!        | USE A STANDARD TEST AS IN THE QR ALGORITHM |
!        | REFERENCE: LAPACK SUBROUTINE ZLAHQR.       |
!        | NOTE: SINCE THE SUBDIAGONALS OF THE        |
!        | COMPRESSED H ARE NONNEGATIVE REAL NUMBERS, |
!        | WE TAKE ADVANTAGE OF THIS.                 |
!        %--------------------------------------------%
!
        tst1 = zabs1(h(i, i))+zabs1(h(i+1, i+1))
        b_lda = to_blas_int(ldh)
        b_n = to_blas_int(kev)
        if (tst1 .eq. rzero) tst1 = zlanhs('1', b_n, h, b_lda, rbid)
        if (dble(h(i+1, i)) .le. max(ulp*tst1, smlnum)) h(i+1, i) = zero
    end do
!
!     %-------------------------------------------------%
!     | COMPUTE THE (KEV+1)-ST COLUMN OF (V*Q) AND      |
!     | TEMPORARILY STORE THE RESULT IN WORKD(N+1:2*N). |
!     | THIS IS NEEDED IN THE RESIDUAL UPDATE SINCE WE  |
!     | CANNOT GUARANTEE THAT THE CORRESPONDING ENTRY   |
!     | OF H WOULD BE ZERO AS IN EXACT ARITHMETIC.      |
!     %-------------------------------------------------%
!
    if (dble(h(kev+1, kev)) .gt. rzero) then
        b_lda = to_blas_int(ldv)
        b_m = to_blas_int(n)
        b_n = to_blas_int(kplusp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zgemv('N', b_m, b_n, one, v, &
                   b_lda, q(1, kev+1), b_incx, zero, workd(n+1), &
                   b_incy)
    end if
!
!     %----------------------------------------------------------%
!     | COMPUTE COLUMN 1 TO KEV OF (V*Q) IN BACKWARD ORDER       |
!     | TAKING ADVANTAGE OF THE UPPER HESSENBERG STRUCTURE OF Q. |
!     %----------------------------------------------------------%
!
    do i = 1, kev
        b_lda = to_blas_int(ldv)
        b_m = to_blas_int(n)
        b_n = to_blas_int(kplusp-i+1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zgemv('N', b_m, b_n, one, v, &
                   b_lda, q(1, kev-i+1), b_incx, zero, workd, &
                   b_incy)
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, workd, b_incx, v(1, kplusp-i+1), b_incy)
    end do
!
!     %-------------------------------------------------%
!     |  MOVE V(:,KPLUSP-KEV+1:KPLUSP) INTO V(:,1:KEV). |
!     %-------------------------------------------------%
!
    b_ldb = to_blas_int(ldv)
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(kev)
    call zlacpy('A', b_m, b_n, v(1, kplusp-kev+1), b_lda, &
                v, b_ldb)
!
!     %--------------------------------------------------------------%
!     | COPY THE (KEV+1)-ST COLUMN OF (V*Q) IN THE APPROPRIATE PLACE |
!     %--------------------------------------------------------------%
!
    if (dble(h(kev+1, kev)) .gt. rzero) then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, workd(n+1), b_incx, v(1, kev+1), b_incy)
    end if
!
!     %-------------------------------------%
!     | UPDATE THE RESIDUAL VECTOR:         |
!     |    R <- SIGMAK*R + BETAK*V(:,KEV+1) |
!     | WHERE                               |
!     |    SIGMAK = (E_{KEV+P}'*Q)*E_{KEV}  |
!     |    BETAK = E_{KEV+1}'*H*E_{KEV}     |
!     %-------------------------------------%
!
    call zlscal(n, q(kplusp, kev), resid, 1)
    if (dble(h(kev+1, kev)) .gt. rzero) then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zaxpy(b_n, h(kev+1, kev), v(1, kev+1), b_incx, resid, &
                   b_incy)
    end if
!
    if (msglvl .gt. 1) then
        call zvout(logfil, 1, q(kplusp, kev), ndigit, '_NAPPS: SIGMAK = (E_(KEV+P)T*Q)*E_(KEV)')
        call zvout(logfil, 1, h(kev+1, kev), ndigit, '_NAPPS: BETAK = E_(KEV+1)T*H*E_(KEV)')
        call ivout(logfil, 1, [kev], ndigit, '_NAPPS: ORDER OF THE FINAL HESSENBERG MATRIX ')
        if (msglvl .gt. 2) then
            call zmout(logfil, kev, kev, h, ldh, &
                       ndigit, '_NAPPS: UPDATED HESSENBERG MATRIX H FOR NEXT ITERATION')
        end if
!
    end if
!
9000 continue
    call matfpe(1)
!
!
!     %---------------%
!     | END OF ZNAPPS |
!     %---------------%
!
end subroutine
