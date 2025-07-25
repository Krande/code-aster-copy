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
subroutine zgetv0(ido, bmat, initv, n, j, &
                  v, ldv, resid, rnorm, ipntr, &
                  workd, ierr, alpha)
!
!     SUBROUTINE ARPACK GENERANT UN VECTEUR INITIAL DANS IM(OP).
!-----------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZGETV0
!
!\DESCRIPTION:
!  GENERATE A RANDOM INITIAL RESIDUAL VECTOR FOR THE ARNOLDI PROCESS.
!  FORCE THE RESIDUAL VECTOR TO BE IN THE RANGE OF THE OPERATOR OP.
!
!\USAGE:
!  CALL ZGETV0
!     ( IDO, BMAT, INITV, N, J, V, LDV, RESID, RNORM,
!       IPNTR, WORKD, IERR )
!
!\ARGUMENTS
!  IDO     INTEGER.  (INPUT/OUTPUT)
!          REVERSE COMMUNICATION FLAG.  IDO MUST BE ZERO ON THE FIRST
!          CALL TO ZGETV0.
!          -------------------------------------------------------------
!          IDO =  0: FIRST CALL TO THE REVERSE COMMUNICATION INTERFACE
!          IDO = -1: COMPUTE  Y = OP * X  WHERE
!                    IPNTR(1) IS THE POINTER INTO WORKD FOR X,
!                    IPNTR(2) IS THE POINTER INTO WORKD FOR Y.
!                    THIS IS FOR THE INITIALIZATION PHASE TO FORCE THE
!                    STARTING VECTOR INTO THE RANGE OF OP.
!          IDO =  2: COMPUTE  Y = B * X  WHERE
!                    IPNTR(1) IS THE POINTER INTO WORKD FOR X,
!                    IPNTR(2) IS THE POINTER INTO WORKD FOR Y.
!          IDO = 99: DONE
!         -------------------------------------------------------------
!
!  BMAT    CHARACTER*1.  (INPUT)
!          BMAT SPECIFIES THE TYPE OF THE MATRIX B IN THE (GENERALIZED)
!          EIGENVALUE PROBLEM A*X = LAMBDA*B*X.
!          B = 'I' -> STANDARD EIGENVALUE PROBLEM A*X = LAMBDA*X
!          B = 'G' -> GENERALIZED EIGENVALUE PROBLEM A*X = LAMBDA*B*X
!
!
!  INITV   LOGICAL VARIABLE.  (INPUT)
!          .TRUE.  => THE INITIAL RESIDUAL VECTOR IS GIVEN IN RESID.
!          .FALSE. => GENERATE A RANDOM INITIAL RESIDUAL VECTOR.
!
!  N       INTEGER.  (INPUT)
!          DIMENSION OF THE PROBLEM.
!
!  J       INTEGER.  (INPUT)
!          INDEX OF THE RESIDUAL VECTOR TO BE GENERATED, WITH RESPECT
!          TO THE ARNOLDI PROCESS.  J > 1 IN CASE OF A "RESTART".
!
!  V       COMPLEX*16 N BY J ARRAY.  (INPUT)
!          THE FIRST J-1 COLUMNS OF V CONTAIN THE CURRENT ARNOLDI BASIS
!          IF THIS IS A "RESTART".
!
!  LDV     INTEGER.  (INPUT)
!          LEADING DIMENSION OF V EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  RESID   COMPLEX*16 ARRAY OF LENGTH N.  (INPUT/OUTPUT)
!          INITIAL RESIDUAL VECTOR TO BE GENERATED.  IF RESID IS
!          PROVIDED, FORCE RESID INTO THE RANGE OF THE OPERATOR OP.
!
!  RNORM   DOUBLE PRECISION SCALAR.  (OUTPUT)
!          B-NORM OF THE GENERATED RESIDUAL.
!
!  IPNTR   INTEGER ARRAY OF LENGTH 3.  (OUTPUT)
!
!  WORKD   COMPLEX*16 WORK ARRAY OF LENGTH 2*N.(REVERSE COMMUNICATION).
!          ON EXIT, WORK(1:N) = B*RESID TO BE USED IN SSAITR.
!
!  IERR    INTEGER.  (OUTPUT)
!          =  0: NORMAL EXIT.
!          = -1: CANNOT GENERATE A NONTRIVIAL RESTARTED RESIDUAL VECTOR
!                IN THE RANGE OF THE OPERATOR OP.
!
!\ENDDOC
!
!----------------------------------------------------------------------
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
!   ZVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!   ZLARNV  LAPACK ROUTINE FOR GENERATING A RANDOM VECTOR.
!   ZGEMV   LEVEL 2 BLAS ROUTINE FOR MATRIX VECTOR MULTIPLICATION.
!   ZCOPY   LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER.
!   ZDOTC   LEVEL 1 BLAS THAT COMPUTES THE SCALAR PRODUCT OF TWO
!            VECTORS.
!   DZNRM2  LEVEL 1 BLAS THAT COMPUTES THE NORM OF A VECTOR.
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
! FILE: GETV0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     %-----------------------------%
!     | INCLUDE FILES FOR DEBUGGING |
!     %-----------------------------%
!
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/dvout.h"
#include "asterfort/zvout.h"
#include "blas/dlapy2.h"
#include "blas/dznrm2.h"
#include "blas/zcopy.h"
#include "blas/zdotc.h"
#include "blas/zgemv.h"
#include "blas/zlarnv.h"
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr, mneigh, mnapps
    integer(kind=8) :: mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
    integer(kind=8) :: nopx, nbx, nrorth, nitref, nrstrt
    common/infor/&
     &  nopx, nbx, nrorth, nitref, nrstrt
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
    character(len=1) :: bmat
    aster_logical :: initv
    integer(kind=8) :: ido, ierr, j, ldv, n
    real(kind=8) :: rnorm, alpha
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    integer(kind=8) :: ipntr(3)
    complex(kind=8) :: resid(n), v(ldv, j), workd(2*n)
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    complex(kind=8) :: one, zero
    real(kind=8) :: rzero
    parameter(one=(1.0d+0, 0.0d+0), zero=(0.0d+0, 0.0d+0),&
     &            rzero=0.0d+0)
!
!     %------------------------%
!     | LOCAL SCALARS & ARRAYS |
!     %------------------------%
!
    aster_logical :: first, inits, orth
    integer(kind=4) :: iseed4(4)
    integer(kind=8) :: idist, iseed(4), iter, msglvl, jj
    real(kind=8) :: rnorm0
    complex(kind=8) :: cnorm
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda, b_m
    blas_int :: b_idist
    save first, iseed, inits, iter, msglvl, orth, rnorm0
!
!     %--------------------%
!     | EXTERNAL FUNCTIONS |
!     %--------------------%
!
!
!     %-----------------%
!     | DATA STATEMENTS |
!     %-----------------%
!
    data inits/.true./
!
!     %-----------------------%
!     | EXECUTABLE STATEMENTS |
!     %-----------------------%
!
    call matfpe(-1)
!
!     %-----------------------------------%
!     | INITIALIZE THE SEED OF THE LAPACK |
!     | RANDOM NUMBER GENERATOR           |
!     %-----------------------------------%
!
    if (inits) then
        iseed(1) = 1
        iseed(2) = 3
        iseed(3) = 5
        iseed(4) = 7
        inits = .false.
    end if
!
    if (ido .eq. 0) then
!
!        %-------------------------------%
!        | INITIALIZE TIMING STATISTICS  |
!        | & MESSAGE LEVEL FOR DEBUGGING |
!        %-------------------------------%
!
        msglvl = mgetv0
        ierr = 0
        iter = 0
        first = .false.
        orth = .false.
!
!        %-----------------------------------------------------%
!        | POSSIBLY GENERATE A RANDOM STARTING VECTOR IN RESID |
!        | USE A LAPACK RANDOM NUMBER GENERATOR USED BY THE    |
!        | MATRIX GENERATION ROUTINES.                         |
!        |    IDIST = 1: UNIFORM (0,1)  DISTRIBUTION;          |
!        |    IDIST = 2: UNIFORM (-1,1) DISTRIBUTION;          |
!        |    IDIST = 3: NORMAL  (0,1)  DISTRIBUTION;          |
!        %-----------------------------------------------------%
!
        if (.not. initv) then
            idist = 2
            iseed4(1) = iseed(1)
            iseed4(2) = iseed(2)
            iseed4(3) = iseed(3)
            iseed4(4) = iseed(4)
            b_idist = to_blas_int(idist)
            b_n = to_blas_int(n)
            call zlarnv(b_idist, iseed4, b_n, resid)
            iseed(1) = iseed4(1)
            iseed(2) = iseed4(2)
            iseed(3) = iseed4(3)
            iseed(4) = iseed4(4)
        end if
!
!        %----------------------------------------------------------%
!        | FORCE THE STARTING VECTOR INTO THE RANGE OF OP TO HANDLE |
!        | THE GENERALIZED PROBLEM WHEN B IS POSSIBLY (SINGULAR).   |
!        %----------------------------------------------------------%
!
        if (bmat .eq. 'G') then
            nopx = nopx+1
            ipntr(1) = 1
            ipntr(2) = n+1
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, resid, b_incx, workd, b_incy)
            ido = -1
            goto 9000
        end if
    end if
!
!     %----------------------------------------%
!     | BACK FROM COMPUTING B*(INITIAL-VECTOR) |
!     %----------------------------------------%
!
    if (first) goto 20
!
!     %-----------------------------------------------%
!     | BACK FROM COMPUTING B*(ORTHOGONALIZED-VECTOR) |
!     %-----------------------------------------------%
!
    if (orth) goto 40
!
!     %------------------------------------------------------%
!     | STARTING VECTOR IS NOW IN THE RANGE OF OP; R = OP*R; |
!     | COMPUTE B-NORM OF STARTING VECTOR.                   |
!     %------------------------------------------------------%
!
    first = .true.
    if (bmat .eq. 'G') then
        nbx = nbx+1
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, workd(n+1), b_incx, resid, b_incy)
        ipntr(1) = n+1
        ipntr(2) = 1
        ido = 2
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd, b_incy)
    end if
!
20  continue
!
    first = .false.
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cnorm = zdotc(b_n, resid, b_incx, workd, b_incy)
        rnorm0 = sqrt(dlapy2(dble(cnorm), dimag(cnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        rnorm0 = dznrm2(b_n, resid, b_incx)
    end if
    rnorm = rnorm0
!
!     %---------------------------------------------%
!     | EXIT IF THIS IS THE VERY FIRST ARNOLDI STEP |
!     %---------------------------------------------%
!
    if (j .eq. 1) goto 50
!
!    %----------------------------------------------------------------
!    | OTHERWISE NEED TO B-ORTHOGONALIZE THE STARTING VECTOR AGAINST |
!    | THE CURRENT ARNOLDI BASIS USING GRAM-SCHMIDT WITH ITER. REF.  |
!    | THIS IS THE CASE WHERE AN INVARIANT SUBSPACE IS ENCOUNTERED   |
!    | IN THE MIDDLE OF THE ARNOLDI FACTORIZATION.                   |
!    |                                                               |
!    |       S = V^{T}*B*R;   R = R - V*S;                           |
!    |                                                               |
!    | STOPPING CRITERIA USED FOR ITER. REF. IS DISCUSSED IN         |
!    | PARLETT'S BOOK, PAGE 107 AND IN GRAGG & REICHEL TOMS PAPER.   |
!    %---------------------------------------------------------------%
!
    orth = .true.
30  continue
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('C', b_m, b_n, one, v, &
               b_lda, workd, b_incx, zero, workd(n+1), &
               b_incy)
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('N', b_m, b_n, -one, v, &
               b_lda, workd(n+1), b_incx, one, resid, &
               b_incy)
!
!     %----------------------------------------------------------%
!     | COMPUTE THE B-NORM OF THE ORTHOGONALIZED STARTING VECTOR |
!     %----------------------------------------------------------%
!
    if (bmat .eq. 'G') then
        nbx = nbx+1
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(n+1), b_incy)
        ipntr(1) = n+1
        ipntr(2) = 1
        ido = 2
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd, b_incy)
    end if
!
40  continue
!
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cnorm = zdotc(b_n, resid, b_incx, workd, b_incy)
        rnorm = sqrt(dlapy2(dble(cnorm), dimag(cnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        rnorm = dznrm2(b_n, resid, b_incx)
    end if
!
!     %--------------------------------------%
!     | CHECK FOR FURTHER ORTHOGONALIZATION. |
!     %--------------------------------------%
!
    if (msglvl .gt. 2) then
        call dvout(logfil, 1, [rnorm0], ndigit, '_GETV0: RE-ORTHONALIZATION ; RNORM0 IS')
        call dvout(logfil, 1, [rnorm], ndigit, '_GETV0: RE-ORTHONALIZATION ; RNORM IS')
    end if
!
    if (rnorm .gt. alpha*rnorm0) goto 50
!
    iter = iter+1
    if (iter .le. 1) then
!
!        %-----------------------------------%
!        | PERFORM ITERATIVE REFINEMENT STEP |
!        %-----------------------------------%
!
        rnorm0 = rnorm
        goto 30
    else
!
!        %------------------------------------%
!        | ITERATIVE REFINEMENT STEP "FAILED" |
!        %------------------------------------%
!
        do jj = 1, n
            resid(jj) = zero
        end do
        rnorm = rzero
        ierr = -1
    end if
!
50  continue
!
    if (msglvl .gt. 0) then
        call dvout(logfil, 1, [rnorm], ndigit, &
                   '_GETV0: B-NORM OF INITIAL / RESTARTED STARTING VECTOR')
    end if
    if (msglvl .gt. 2) then
        call zvout(logfil, n, resid, ndigit, '_GETV0: INITIAL / RESTARTED STARTING VECTOR')
    end if
    ido = 99
!
9000 continue
    call matfpe(1)
!
!     %---------------%
!     | END OF ZGETV0 |
!     %---------------%
!
end subroutine
