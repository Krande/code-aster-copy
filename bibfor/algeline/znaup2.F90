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
subroutine znaup2(ido, bmat, n, which, nev, &
                  np, tol, resid, ishift, mxiter, &
                  v, ldv, h, ldh, ritz, &
                  bounds, q, ldq, workl, ipntr, &
                  workd, rwork, info, neqact, alpha)
!
!     SUBROUTINE ARPACK CALCULANT LES VALEURS PROPRES DE (OP) VIA
!     IRAM.
!---------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZNAUP2
!
!\DESCRIPTION:
!  INTERMEDIATE LEVEL INTERFACE CALLED BY ZNAUPD .
!
!\USAGE:
!  CALL ZNAUP2
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS,
!       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
!
!\ARGUMENTS
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: SAME AS DEFINED IN ZNAUPD .
!  ISHIFT, MXITER: SEE THE DEFINITION OF IPARAM IN ZNAUPD .
!
!  NP      INTEGER.  (INPUT/OUTPUT)
!          CONTAINS THE NUMBER OF IMPLICIT SHIFTS TO APPLY DURING
!          EACH ARNOLDI ITERATION.
!          IF ISHIFT=1, NP IS ADJUSTED DYNAMICALLY AT EACH ITERATION
!          TO ACCELERATE CONVERGENCE AND PREVENT STAGNATION.
!          THIS IS ALSO ROUGHLY EQUAL TO THE NUMBER OF MATRIX-VECTOR
!          PRODUCTS (INVOLVING THE OPERATOR OP) PER ARNOLDI ITERATION.
!          THE LOGIC FOR ADJUSTING IS CONTAINED WITHIN THE CURRENT
!          SUBROUTINE.
!          IF ISHIFT=0, NP IS THE NUMBER OF SHIFTS THE USER NEEDS
!          TO PROVIDE VIA REVERSE COMUNICATION. 0 < NP < NCV-NEV.
!          NP MAY BE LESS THAN NCV-NEV SINCE A LEADING BLOCK OF THE
!          CURRENT
!          UPPER HESSENBERG MATRIX HAS SPLIT OFF AND CONTAINS "UNWANTED"
!          RITZ VALUES.
!          UPON TERMINATION OF THE IRA ITERATION, NP CONTAINS THE NUMBER
!          OF "CONVERGED" WANTED RITZ VALUES.
!
!  V       COMPLEX*16  N BY (NEV+NP) ARRAY.  (INPUT/OUTPUT)
!          THE ARNOLDI BASIS VECTORS ARE RETURNED IN THE FIRST NEV
!          COLUMNS OF V.
!
!  LDV     INTEGER.  (INPUT)
!          LEADING DIMENSION OF V EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  H       COMPLEX*16  (NEV+NP) BY (NEV+NP) ARRAY.  (OUTPUT)
!          H IS USED TO STORE THE GENERATED UPPER HESSENBERG MATRIX
!
!  LDH     INTEGER.  (INPUT)
!          LEADING DIMENSION OF H EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  RITZ    COMPLEX*16  ARRAY OF LENGTH NEV+NP.  (OUTPUT)
!          RITZ(1:NEV)  CONTAINS THE COMPUTED RITZ VALUES OF OP.
!
!  BOUNDS  COMPLEX*16  ARRAY OF LENGTH NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) CONTAIN THE ERROR BOUNDS CORRESPONDING TO
!          THE COMPUTED RITZ VALUES.
!
!  Q       COMPLEX*16  (NEV+NP) BY (NEV+NP) ARRAY.  (WORKSPACE)
!          PRIVATE (REPLICATED) WORK ARRAY USED TO ACCUMULATE THE
!          ROTATION IN THE SHIFT APPLICATION STEP.
!
!  LDQ     INTEGER.  (INPUT)
!          LEADING DIMENSION OF Q EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  WORKL   COMPLEX*16  WORK ARRAY OF LENGTH AT LEAST
!          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.  IT IS USED IN SHIFTS CALCULATION, SHIFTS
!          APPLICATION AND CONVERGENCE CHECKING.
!
!
!  IPNTR   INTEGER ARRAY OF LENGTH 3.  (OUTPUT)
!          POINTER TO MARK THE STARTING LOCATIONS IN THE WORKD FOR
!          VECTORS USED BY THE ARNOLDI ITERATION.
!          -------------------------------------------------------------
!          IPNTR(1): POINTER TO THE CURRENT OPERAND VECTOR X.
!          IPNTR(2): POINTER TO THE CURRENT RESULT VECTOR Y.
!          IPNTR(3): POINTER TO THE VECTOR B * X WHEN USED IN THE
!                    SHIFT-AND-INVERT MODE.  X IS THE CURRENT OPERAND.
!          -------------------------------------------------------------
!
!  WORKD   COMPLEX*16  WORK ARRAY OF LENGTH 3*N.  (WORKSPACE)
!          DISTRIBUTED ARRAY TO BE USED IN THE BASIC ARNOLDI ITERATION
!          FOR REVERSE COMMUNICATION.  THE USER SHOULD NOT USE WORKD
!          AS TEMPORARY WORKSPACE DURING THE ITERATION !!!!!!!!!!
!          SEE DATA DISTRIBUTION NOTE IN ZNAUPD .
!
!  RWORK   DOUBLE PRECISION    WORK ARRAY OF LENGTH  NEV+NP ( WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.
!
!  INFO    INTEGER.  (INPUT/OUTPUT)
!          IF INFO .EQ. 0, A RANDOMLY INITIAL RESIDUAL VECTOR IS USED.
!          IF INFO .NE. 0, RESID CONTAINS THE INITIAL RESIDUAL VECTOR,
!                          POSSIBLY FROM A PREVIOUS RUN.
!          ERROR FLAG ON OUTPUT.
!          =     0: NORMAL RETURN.
!          =     1: MAXIMUM NUMBER OF ITERATIONS TAKEN.
!                   ALL POSSIBLE EIGENVALUES OF OP HAS BEEN FOUND.
!                   NP RETURNS THE NUMBER OF CONVERGED RITZ VALUES.
!          =     2: NO SHIFTS COULD BE APPLIED.
!          =    -8: ERROR RETURN FROM LAPACK EIGENVALUE CALCULATION;
!                   THIS SHOULD NEVER HAPPEN.
!          =    -9: STARTING VECTOR IS ZERO.
!          = -9999: COULD NOT BUILD AN ARNOLDI FACTORIZATION.
!                   SIZE THAT WAS BUILT IN RETURNED IN NP.
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
!  2. R.B. LEHOUCQ, "ANALYSIS AND IMPLEMENTATION OF AN IMPLICITLY
!     RESTARTED ARNOLDI ITERATION", RICE UNIVERSITY TECHNICAL REPORT
!     TR95-13, DEPARTMENT OF COMPUTATIONAL AND APPLIED MATHEMATICS.
!
!\ROUTINES CALLED:
!     ZGETV0   ARPACK INITIAL VECTOR GENERATION ROUTINE.
!     ZNAITR   ARPACK ARNOLDI FACTORIZATION ROUTINE.
!     ZNAPPS   ARPACK APPLICATION OF IMPLICIT SHIFTS ROUTINE.
!     ZNEIGH   ARPACK COMPUTE RITZ VALUES AND ERROR BOUNDS ROUTINE.
!     ZNGETS   ARPACK REORDER RITZ VALUES AND ERROR BOUNDS ROUTINE.
!     ZSORTC   ARPACK SORTING ROUTINE.
!     IVOUT   ARPACK UTILITY ROUTINE THAT PRINTS INTEGERS.
!     ZMOUT    ARPACK UTILITY ROUTINE THAT PRINTS MATRICES
!     ZVOUT    ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     DVOUT    ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     DLAPY2   LAPACK ROUTINE TO COMPUTE SQRT(X**2+Y**2) CAREFULLY.
!     ZCOPY    LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER .
!     ZDOTC    LEVEL 1 BLAS THAT COMPUTES THE SCALAR PRODUCT OF TWO
!               VECTORS.
!     GLSWAP    LEVEL 1 BLAS THAT SWAPS TWO VECTORS.
!     DZNRM2   LEVEL 1 BLAS THAT COMPUTES THE NORM OF A VECTOR.
!
!\AUTHOR
!     DANNY SORENSEN               PHUONG VU
!     RICHARD LEHOUCQ              CRPC / RICE UNIVERSITYA
!     CHAO YANG                    HOUSTON, TEXAS
!     DEPT. OF COMPUTATIONAL&
!     APPLIED MATHEMATICS
!     RICE UNIVERSITY
!     HOUSTON, TEXAS
!
!\SCCS INFORMATION: @(#)
! FILE: NAUP2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
!
!\REMARKS
!     1. NONE
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1504
    implicit none
!
!     %----------------------------------------------------%
!     | INCLUDE FILES FOR DEBUGGING AND TIMING INFORMATION |
!     %----------------------------------------------------%
!
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterc/r8prem.h"
#include "asterfort/dvout.h"
#include "asterfort/ivout.h"
#include "asterfort/zgetv0.h"
#include "asterfort/zmout.h"
#include "asterfort/znaitr.h"
#include "asterfort/znapps.h"
#include "asterfort/zneigh.h"
#include "asterfort/zngets.h"
#include "asterfort/zsortc.h"
#include "asterfort/zvout.h"
#include "blas/dlapy2.h"
#include "blas/dznrm2.h"
#include "blas/zcopy.h"
#include "blas/zdotc.h"
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr, mneigh, mnapps
    integer(kind=8) :: mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
    integer(kind=8) :: nopx, nbx, nrorth, nitref, nrstrt
    common/infor/&
     &  nopx, nbx, nrorth, nitref, nrstrt
!
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
    character(len=1) :: bmat
    character(len=2) :: which
    integer(kind=8) :: ido, info, ishift, ldh, ldq, ldv, mxiter, n, nev, np, neqact
    real(kind=8) :: tol, alpha
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    integer(kind=8) :: ipntr(13)
    complex(kind=8) :: bounds(nev+np), h(ldh, nev+np), q(ldq, nev+np), resid(n)
    complex(kind=8) :: ritz(nev+np), v(ldv, nev+np), workd(3*n)
    complex(kind=8) :: workl((nev+np)*(nev+np+3))
    real(kind=8) :: rwork(nev+np)
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    complex(kind=8) :: zero
    real(kind=8) :: rzero
    parameter(zero=(0.0d+0, 0.0d+0), rzero=0.0d+0)
!
!     %---------------%
!     | LOCAL SCALARS |
!     %---------------%
!
    aster_logical :: cnorm, getv0, initv, update, ushift
    integer(kind=8) :: ierr, iter, kplusp, msglvl, nconv, nevbef, nev0, np0, nptemp, i
    integer(kind=8) :: j
    complex(kind=8) :: cpnorm
    real(kind=8) :: rnorm, eps23, rtemp
    character(len=2) :: wprime
!
    save cnorm, getv0, initv, update, ushift,&
     &           rnorm, iter, kplusp, msglvl, nconv,&
     &           nevbef, nev0, np0, eps23
!
!
!     %-----------------------%
!     | LOCAL ARRAY ARGUMENTS |
!     %-----------------------%
!
    integer(kind=8) :: kp(3)
    blas_int :: b_incx, b_incy, b_n
!
!     %--------------------%
!     | EXTERNAL FUNCTIONS |
!     %--------------------%
!
!
!     %-----------------------%
!     | EXECUTABLE STATEMENTS |
!     %-----------------------%
!
    call matfpe(-1)
!
    if (ido .eq. 0) then
!
        msglvl = mnaup2
!        %-------------------------------------%
!        | GET THE MACHINE DEPENDENT CONSTANT. |
!        %-------------------------------------%
!
        eps23 = r8prem()*0.5d0
        eps23 = eps23**(2.0d+0/3.0d+0)
!
        nev0 = nev
        np0 = np
!
!        %-------------------------------------%
!        | KPLUSP IS THE BOUND ON THE LARGEST  |
!        |        LANCZOS FACTORIZATION BUILT. |
!        | NCONV IS THE CURRENT NUMBER OF      |
!        |        "CONVERGED" EIGENVALUES.     |
!        | ITER IS THE COUNTER ON THE CURRENT  |
!        |      ITERATION STEP.                |
!        %-------------------------------------%
!
        kplusp = nev+np
        nconv = 0
        iter = 0
!
!        %---------------------------------%
!        | GET MACHINE DEPENDENT CONSTANT. |
!        %---------------------------------%
!
        eps23 = r8prem()*0.5d0
        eps23 = eps23**(2.0d+0/3.0d+0)
!
!        %---------------------------------------%
!        | SET FLAGS FOR COMPUTING THE FIRST NEV |
!        | STEPS OF THE ARNOLDI FACTORIZATION.   |
!        %---------------------------------------%
!
        getv0 = .true.
        update = .false.
        ushift = .false.
        cnorm = .false.
!
        if (info .ne. 0) then
!
!           %--------------------------------------------%
!           | USER PROVIDES THE INITIAL RESIDUAL VECTOR. |
!           %--------------------------------------------%
!
            initv = .true.
            info = 0
        else
            initv = .false.
        end if
    end if
!
!     %---------------------------------------------%
!     | GET A POSSIBLY RANDOM STARTING VECTOR AND   |
!     | FORCE IT INTO THE RANGE OF THE OPERATOR OP. |
!     %---------------------------------------------%
!
!
    if (getv0) then
        call zgetv0(ido, bmat, initv, n, 1, &
                    v, ldv, resid, rnorm, ipntr, &
                    workd, info, alpha)
!
        if (ido .ne. 99) goto 9000
!
        if (rnorm .eq. rzero) then
!
!           %-----------------------------------------%
!           | THE INITIAL VECTOR IS ZERO. ERROR EXIT. |
!           %-----------------------------------------%
!
            info = -9
            goto 1100
        end if
        getv0 = .false.
        ido = 0
    end if
!
!     %-----------------------------------%
!     | BACK FROM REVERSE COMMUNICATION : |
!     | CONTINUE WITH UPDATE STEP         |
!     %-----------------------------------%
!
    if (update) goto 20
!
!     %-------------------------------------------%
!     | BACK FROM COMPUTING USER SPECIFIED SHIFTS |
!     %-------------------------------------------%
!
    if (ushift) goto 50
!
!     %-------------------------------------%
!     | BACK FROM COMPUTING RESIDUAL NORM   |
!     | AT THE END OF THE CURRENT ITERATION |
!     %-------------------------------------%
!
    if (cnorm) goto 100
!
!     %----------------------------------------------------------%
!     | COMPUTE THE FIRST NEV STEPS OF THE ARNOLDI FACTORIZATION |
!     %----------------------------------------------------------%
!
    call znaitr(ido, bmat, n, 0, nev, &
                resid, rnorm, v, ldv, h, &
                ldh, ipntr, workd, info, alpha)
!
    if (ido .ne. 99) goto 9000
!
    if (info .gt. 0) then
        np = info
        mxiter = iter
        info = -9999
        goto 1200
    end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
!     |           EACH ITERATION IMPLICITLY RESTARTS THE ARNOLDI     |
!     |           FACTORIZATION IN PLACE.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
!
1000 continue
!
    iter = iter+1
!
    if (msglvl .gt. 0) then
        call ivout(logfil, 1, [iter], ndigit, &
                   '_NAUP2: **** START OF MAJOR ITERATION NUMBER ****')
    end if
!
!        %-----------------------------------------------------------%
!        | COMPUTE NP ADDITIONAL STEPS OF THE ARNOLDI FACTORIZATION. |
!        | ADJUST NP SINCE NEV MIGHT HAVE BEEN UPDATED BY LAST CALL  |
!        | TO THE SHIFT APPLICATION ROUTINE ZNAPPS .                  |
!        %-----------------------------------------------------------%
!
    np = kplusp-nev
!
    if (msglvl .gt. 1) then
        call ivout(logfil, 1, [nev], ndigit, &
                   '_NAUP2: THE LENGTH OF THE CURRENT ARNOLDI FACTORIZATION')
        call ivout(logfil, 1, [np], ndigit, '_NAUP2: EXTEND THE ARNOLDI FACTORIZATION BY')
    end if
!
!        %-----------------------------------------------------------%
!        | COMPUTE NP ADDITIONAL STEPS OF THE ARNOLDI FACTORIZATION. |
!        %-----------------------------------------------------------%
!
    ido = 0
20  continue
    update = .true.
!
    call znaitr(ido, bmat, n, nev, np, &
                resid, rnorm, v, ldv, h, &
                ldh, ipntr, workd, info, alpha)
!
    if (ido .ne. 99) goto 9000
!
    if (info .gt. 0) then
        np = info
        mxiter = iter
        if (info .ge. neqact) then
            if (msglvl .gt. 0) then
                write (logfil, *)
                write (logfil, *) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
                write (logfil, *) '& ESPACE INVARIANT DE TAILLE &'
                write (logfil, *) '& NEQACT = ', neqact
                write (logfil, *) '& SHUNTAGE PARTIEL DE ZNAUP2 &'
                write (logfil, *) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
                write (logfil, *)
            end if
        else
            info = -9999
            goto 1200
        end if
        goto 1200
    end if
    update = .false.
!
    if (msglvl .gt. 1) then
        call dvout(logfil, 1, [rnorm], ndigit, '_NAUP2: CORRESPONDING B-NORM OF THE RESIDUAL')
    end if
!
!        %--------------------------------------------------------%
!        | COMPUTE THE EIGENVALUES AND CORRESPONDING ERROR BOUNDS |
!        | OF THE CURRENT UPPER HESSENBERG MATRIX.                |
!        %--------------------------------------------------------%
!
    call zneigh(rnorm, kplusp, h, ldh, ritz, &
                bounds, q, ldq, workl, rwork, &
                ierr)
!
    if (ierr .ne. 0) then
        info = -8
        goto 1200
    end if
!
!        %---------------------------------------------------%
!        | SELECT THE WANTED RITZ VALUES AND THEIR BOUNDS    |
!        | TO BE USED IN THE CONVERGENCE TEST.               |
!        | THE WANTED PART OF THE SPECTRUM AND CORRESPONDING |
!        | ERROR BOUNDS ARE IN THE LAST NEV LOC. OF RITZ,    |
!        | AND BOUNDS RESPECTIVELY.                          |
!        %---------------------------------------------------%
!
    nev = nev0
    np = np0
!
!        %--------------------------------------------------%
!        | MAKE A COPY OF RITZ VALUES AND THE CORRESPONDING |
!        | RITZ ESTIMATES OBTAINED FROM ZNEIGH .             |
!        %--------------------------------------------------%
!
    b_n = to_blas_int(kplusp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, ritz, b_incx, workl(kplusp**2+1), b_incy)
    b_n = to_blas_int(kplusp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, bounds, b_incx, workl(kplusp**2+kplusp+1), b_incy)
!
!        %---------------------------------------------------%
!        | SELECT THE WANTED RITZ VALUES AND THEIR BOUNDS    |
!        | TO BE USED IN THE CONVERGENCE TEST.               |
!        | THE WANTED PART OF THE SPECTRUM AND CORRESPONDING |
!        | BOUNDS ARE IN THE LAST NEV LOC. OF RITZ           |
!        | BOUNDS RESPECTIVELY.                              |
!        %---------------------------------------------------%
!
    call zngets(ishift, which, nev, np, ritz, &
                bounds)
!
!        %------------------------------------------------------------%
!        | CONVERGENCE TEST: CURRENTLY WE USE THE FOLLOWING CRITERIA. |
!        | THE RELATIVE ACCURACY OF A RITZ VALUE IS CONSIDERED        |
!        | ACCEPTABLE IF:                                             |
!        |                                                            |
!        | ERROR_BOUNDS(I) .LE. TOL*MAX(EPS23, MAGNITUDE_OF_RITZ(I)). |
!        |                                                            |
!        %------------------------------------------------------------%
!
    nconv = 0
!
    do i = 1, nev
        rtemp = max(eps23, dlapy2(dble(ritz(np+i)), dimag(ritz(np+i))))
        if (dlapy2(dble(bounds(np+i)), dimag(bounds(np+i))) .le. tol*rtemp) then
            nconv = nconv+1
        end if
    end do
!
    if (msglvl .gt. 2) then
        kp(1) = nev
        kp(2) = np
        kp(3) = nconv
        call ivout(logfil, 3, kp, ndigit, '_NAUP2: NEV, NP, NCONV ARE')
        call zvout(logfil, kplusp, ritz, ndigit, '_NAUP2: THE EIGENVALUES OF H')
        call zvout(logfil, kplusp, bounds, ndigit, &
                   '_NAUP2: RITZ ESTIMATES OF THE CURRENT NCV RITZ VALUES')
    end if
!
!        %---------------------------------------------------------%
!        | COUNT THE NUMBER OF UNWANTED RITZ VALUES THAT HAVE ZERO |
!        | RITZ ESTIMATES. IF ANY RITZ ESTIMATES ARE EQUAL TO ZERO |
!        | THEN A LEADING BLOCK OF H OF ORDER EQUAL TO AT LEAST    |
!        | THE NUMBER OF RITZ VALUES WITH ZERO RITZ ESTIMATES HAS  |
!        | SPLIT OFF. NONE OF THESE RITZ VALUES MAY BE REMOVED BY  |
!        | SHIFTING. DECREASE NP THE NUMBER OF SHIFTS TO APPLY. IF |
!        | NO SHIFTS MAY BE APPLIED, THEN PREPARE TO EXIT          |
!        %---------------------------------------------------------%
!
    nptemp = np
    do j = 1, nptemp
        if (bounds(j) .eq. zero) then
            np = np-1
            nev = nev+1
        end if
    end do
!
    if ((nconv .ge. nev0) .or. (iter .gt. mxiter) .or. (np .eq. 0)) then
!
        if (msglvl .gt. 4) then
            call zvout(logfil, kplusp, workl(kplusp**2+1), ndigit, &
                       '_NAUP2: EIGENVALUES COMPUTED BY _NEIGH:')
            call zvout(logfil, kplusp, workl(kplusp**2+kplusp+1), ndigit, &
                       '_NAUP2: RITZ ESTIMATES COMPUTED BY _NEIGH:')
        end if
!
!           %------------------------------------------------%
!           | PREPARE TO EXIT. PUT THE CONVERGED RITZ VALUES |
!           | AND CORRESPONDING BOUNDS IN RITZ(1:NCONV) AND  |
!           | BOUNDS(1:NCONV) RESPECTIVELY. THEN SORT. BE    |
!           | CAREFUL WHEN NCONV > NP                        |
!           %------------------------------------------------%
!
!           %------------------------------------------%
!           |  USE H( 3,1 ) AS STORAGE TO COMMUNICATE  |
!           |  RNORM TO ZNEUPD  IF NEEDED               |
!           %------------------------------------------%
!
        h(3, 1) = dcmplx(rnorm, rzero)
!
!           %----------------------------------------------%
!           | SORT RITZ VALUES SO THAT CONVERGED RITZ      |
!           | VALUES APPEAR WITHIN THE FIRST NEV LOCATIONS |
!           | OF RITZ AND BOUNDS, AND THE MOST DESIRED ONE |
!           | APPEARS AT THE FRONT.                        |
!           %----------------------------------------------%
!
        if (which .eq. 'LM') wprime = 'SM'
        if (which .eq. 'SM') wprime = 'LM'
        if (which .eq. 'LR') wprime = 'SR'
        if (which .eq. 'SR') wprime = 'LR'
        if (which .eq. 'LI') wprime = 'SI'
        if (which .eq. 'SI') wprime = 'LI'
!
        call zsortc(wprime, .true._1, kplusp, ritz, bounds)
!
!           %--------------------------------------------------%
!           | SCALE THE RITZ ESTIMATE OF EACH RITZ VALUE       |
!           | BY 1 / MAX(EPS23, MAGNITUDE OF THE RITZ VALUE).  |
!           %--------------------------------------------------%
!
        do j = 1, nev0
            rtemp = max(eps23, dlapy2(dble(ritz(j)), dimag(ritz(j))))
            bounds(j) = bounds(j)/rtemp
        end do
!
!           %---------------------------------------------------%
!           | SORT THE RITZ VALUES ACCORDING TO THE SCALED RITZ |
!           | ESTIMATES.  THIS WILL PUSH ALL THE CONVERGED ONES |
!           | TOWARDS THE FRONT OF RITZ, BOUNDS (IN THE CASE    |
!           | WHEN NCONV < NEV.)                                |
!           %---------------------------------------------------%
!
        wprime = 'LM'
        call zsortc(wprime, .true._1, nev0, bounds, ritz)
!
!           %----------------------------------------------%
!           | SCALE THE RITZ ESTIMATE BACK TO ITS ORIGINAL |
!           | VALUE.                                       |
!           %----------------------------------------------%
!
        do j = 1, nev0
            rtemp = max(eps23, dlapy2(dble(ritz(j)), dimag(ritz(j))))
            bounds(j) = bounds(j)*rtemp
        end do
!
!           %-----------------------------------------------%
!           | SORT THE CONVERGED RITZ VALUES AGAIN SO THAT  |
!           | THE "THRESHOLD" VALUE APPEARS AT THE FRONT OF |
!           | RITZ AND BOUND.                               |
!           %-----------------------------------------------%
!
        call zsortc(which, .true._1, nconv, ritz, bounds)
!
        if (msglvl .gt. 1) then
            call zvout(logfil, kplusp, ritz, ndigit, '_NAUP2: SORTED EIGENVALUES')
            call zvout(logfil, kplusp, bounds, ndigit, '_NAUP2: SORTED RITZ ESTIMATES.')
        end if
!
!           %------------------------------------%
!           | MAX ITERATIONS HAVE BEEN EXCEEDED. |
!           %------------------------------------%
!
        if (iter .gt. mxiter .and. nconv .lt. nev0) info = 1
!
!           %---------------------%
!           | NO SHIFTS TO APPLY. |
!           %---------------------%
!
        if (np .eq. 0 .and. nconv .lt. nev0) info = 2
!
        np = nconv
        goto 1100
!
    else if ((nconv .lt. nev0) .and. (ishift .eq. 1)) then
!
!           %-------------------------------------------------%
!           | DO NOT HAVE ALL THE REQUESTED EIGENVALUES YET.  |
!           | TO PREVENT POSSIBLE STAGNATION, ADJUST THE SIZE |
!           | OF NEV.                                         |
!           %-------------------------------------------------%
!
        nevbef = nev
        nev = nev+min(nconv, np/2)
        if (nev .eq. 1 .and. kplusp .ge. 6) then
            nev = kplusp/2
        else if (nev .eq. 1 .and. kplusp .gt. 3) then
            nev = 2
        end if
        np = kplusp-nev
!
!           %---------------------------------------%
!           | IF THE SIZE OF NEV WAS JUST INCREASED |
!           | RESORT THE EIGENVALUES.               |
!           %---------------------------------------%
!
        if (nevbef .lt. nev) call zngets(ishift, which, nev, np, ritz, &
                                         bounds)
!
    end if
!
    if (msglvl .gt. 0) then
        call ivout(logfil, 1, [nconv], ndigit, &
                   '_NAUP2: NO. OF "CONVERGED" RITZ VALUES AT THIS ITER.')
        if (msglvl .gt. 1) then
            kp(1) = nev
            kp(2) = np
            call ivout(logfil, 2, kp, ndigit, '_NAUP2: NEV AND NP ARE')
            call zvout(logfil, nev, ritz(np+1), ndigit, '_NAUP2: "WANTED" RITZ VALUES ')
            call zvout(logfil, nev, bounds(np+1), ndigit, &
                       '_NAUP2: RITZ ESTIMATES OF THE "WANTED" VALUES ')
        end if
    end if
!
    if (ishift .eq. 0) then
!
!           %-------------------------------------------------------%
!           | USER SPECIFIED SHIFTS: POP BACK OUT TO GET THE SHIFTS |
!           | AND RETURN THEM IN THE FIRST 2*NP LOCATIONS OF WORKL. |
!           %-------------------------------------------------------%
!
        ushift = .true.
        ido = 3
        goto 9000
    end if
50  continue
    ushift = .false.
!
    if (ishift .ne. 1) then
!
!            %----------------------------------%
!            | MOVE THE NP SHIFTS FROM WORKL TO |
!            | RITZ, TO FREE UP WORKL           |
!            | FOR NON-EXACT SHIFT CASE.        |
!            %----------------------------------%
!
        b_n = to_blas_int(np)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, workl, b_incx, ritz, b_incy)
    end if
!
    if (msglvl .gt. 2) then
        call ivout(logfil, 1, [np], ndigit, '_NAUP2: THE NUMBER OF SHIFTS TO APPLY ')
        call zvout(logfil, np, ritz, ndigit, '_NAUP2: VALUES OF THE SHIFTS')
        if (ishift .eq. 1) call zvout(logfil, np, bounds, ndigit, &
                                      '_NAUP2: RITZ ESTIMATES OF THE SHIFTS')
    end if
!
!        %---------------------------------------------------------%
!        | APPLY THE NP IMPLICIT SHIFTS BY QR BULGE CHASING.       |
!        | EACH SHIFT IS APPLIED TO THE WHOLE UPPER HESSENBERG     |
!        | MATRIX H.                                               |
!        | THE FIRST 2*N LOCATIONS OF WORKD ARE USED AS WORKSPACE. |
!        %---------------------------------------------------------%
!
    call znapps(n, nev, np, ritz, v, &
                ldv, h, ldh, resid, q, &
                ldq, workl, workd)
!
!        %---------------------------------------------%
!        | COMPUTE THE B-NORM OF THE UPDATED RESIDUAL. |
!        | KEEP B*RESID IN WORKD(1:N) TO BE USED IN    |
!        | THE FIRST STEP OF THE NEXT CALL TO ZNAITR .  |
!        %---------------------------------------------%
!
    cnorm = .true.
    if (bmat .eq. 'G') then
        nbx = nbx+1
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(n+1), b_incy)
        ipntr(1) = n+1
        ipntr(2) = 1
        ido = 2
!
!           %----------------------------------%
!           | EXIT IN ORDER TO COMPUTE B*RESID |
!           %----------------------------------%
!
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd, b_incy)
    end if
!
100 continue
!
!        %----------------------------------%
!        | BACK FROM REVERSE COMMUNICATION; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
!
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cpnorm = zdotc(b_n, resid, b_incx, workd, b_incy)
        rnorm = sqrt(dlapy2(dble(cpnorm), dimag(cpnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        rnorm = dznrm2(b_n, resid, b_incx)
    end if
    cnorm = .false.
!
    if (msglvl .gt. 2) then
        call dvout(logfil, 1, [rnorm], ndigit, &
                   '_NAUP2: B-NORM OF RESIDUAL FOR COMPRESSED FACTORIZATION')
        call zmout(logfil, nev, nev, h, ldh, &
                   ndigit, '_NAUP2: COMPRESSED UPPER HESSENBERG MATRIX H')
    end if
!
    goto 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
1100 continue
!
    mxiter = iter
    nev = nconv
!
1200 continue
    ido = 99
!
9000 continue
    call matfpe(1)
!
!     %---------------%
!     | END OF ZNAUP2  |
!     %---------------%
!
end subroutine
