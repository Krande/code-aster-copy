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
!
!     SUBROUTINE ARPACK OPERANT NP ETAPE D'ARNOLDI A PARTIR D'UNE
!     FACTORISATION D'ORDRE K.
!-----------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZNAITR
!
!\DESCRIPTION:
!  REVERSE COMMUNICATION INTERFACE FOR APPLYING NP ADDITIONAL STEPS TO
!  A K STEP NONSYMMETRIC ARNOLDI FACTORIZATION.
!
!  INPUT:  OP*V_{K}  -  V_{K}*H = R_{K}*E_{K}^T
!
!          WITH (V_{K}^T)*B*V_{K} = I, (V_{K}^T)*B*R_{K} = 0.
!
!  OUTPUT: OP*V_{K+P}  -  V_{K+P}*H = R_{K+P}*E_{K+P}^T
!
!          WITH (V_{K+P}^T)*B*V_{K+P} = I, (V_{K+P}^T)*B*R_{K+P} = 0.
!
!  WHERE OP AND B ARE AS IN ZNAUPD.  THE B-NORM OF R_{K+P} IS ALSO
!  COMPUTED AND RETURNED.
!
!\USAGE:
!  CALL ZNAITR
!     ( IDO, BMAT, N, K, NP, RESID, RNORM, V, LDV, H, LDH,
!       IPNTR, WORKD, INFO )
!
!\ARGUMENTS
!  IDO     INTEGER.  (INPUT/OUTPUT)
!          REVERSE COMMUNICATION FLAG.
!          -------------------------------------------------------------
!          IDO =  0: FIRST CALL TO THE REVERSE COMMUNICATION INTERFACE
!          IDO = -1: COMPUTE  Y = OP * X  WHERE
!                    IPNTR(1) IS THE POINTER INTO WORK FOR X,
!                    IPNTR(2) IS THE POINTER INTO WORK FOR Y.
!                    THIS IS FOR THE RESTART PHASE TO FORCE THE NEW
!                    STARTING VECTOR INTO THE RANGE OF OP.
!          IDO =  1: COMPUTE  Y = OP * X  WHERE
!                    IPNTR(1) IS THE POINTER INTO WORK FOR X,
!                    IPNTR(2) IS THE POINTER INTO WORK FOR Y,
!                    IPNTR(3) IS THE POINTER INTO WORK FOR B * X.
!          IDO =  2: COMPUTE  Y = B * X  WHERE
!                    IPNTR(1) IS THE POINTER INTO WORK FOR X,
!                    IPNTR(2) IS THE POINTER INTO WORK FOR Y.
!          IDO = 99: DONE
!          -------------------------------------------------------------
!          WHEN THE ROUTINE IS USED IN THE "SHIFT-AND-INVERT" MODE, THE
!          VECTOR B * Q IS ALREADY AVAILABLE AND DO NOT NEED TO BE
!          RECOMPUTED IN FORMING OP * Q.
!
!  BMAT    CHARACTER*1.  (INPUT)
!          BMAT SPECIFIES THE TYPE OF THE MATRIX B THAT DEFINES THE
!          SEMI-INNER PRODUCT FOR THE OPERATOR OP.  SEE ZNAUPD.
!          B = 'I' -> STANDARD EIGENVALUE PROBLEM A*X = LAMBDA*X
!          B = 'G' -> GENERALIZED EIGENVALUE PROBLEM A*X = LAMBDA*M**X
!
!  N       INTEGER.  (INPUT)
!          DIMENSION OF THE EIGENPROBLEM.
!
!  K       INTEGER.  (INPUT)
!          CURRENT SIZE OF V AND H.
!
!  NP      INTEGER.  (INPUT)
!          NUMBER OF ADDITIONAL ARNOLDI STEPS TO TAKE.
!
!  RESID   COMPLEX*16 ARRAY OF LENGTH N.  (INPUT/OUTPUT)
!          ON INPUT:  RESID CONTAINS THE RESIDUAL VECTOR R_{K}.
!          ON OUTPUT: RESID CONTAINS THE RESIDUAL VECTOR R_{K+P}.
!
!  RNORM   DOUBLE PRECISION SCALAR.  (INPUT/OUTPUT)
!          B-NORM OF THE STARTING RESIDUAL ON INPUT.
!          B-NORM OF THE UPDATED RESIDUAL R_{K+P} ON OUTPUT.
!
!  V       COMPLEX*16 N BY K+NP ARRAY.  (INPUT/OUTPUT)
!          ON INPUT:  V CONTAINS THE ARNOLDI VECTORS IN THE FIRST K
!          COLUMNS.
!          ON OUTPUT: V CONTAINS THE NEW NP ARNOLDI VECTORS IN THE NEXT
!          NP COLUMNS.  THE FIRST K COLUMNS ARE UNCHANGED.
!
!  LDV     INTEGER.  (INPUT)
!          LEADING DIMENSION OF V EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  H       COMPLEX*16 (K+NP) BY (K+NP) ARRAY.  (INPUT/OUTPUT)
!          H IS USED TO STORE THE GENERATED UPPER HESSENBERG MATRIX.
!
!  LDH     INTEGER.  (INPUT)
!          LEADING DIMENSION OF H EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  IPNTR   INTEGER ARRAY OF LENGTH 3.  (OUTPUT)
!          POINTER TO MARK THE STARTING LOCATIONS IN THE WORK FOR
!          VECTORS USED BY THE ARNOLDI ITERATION.
!          -------------------------------------------------------------
!          IPNTR(1): POINTER TO THE CURRENT OPERAND VECTOR X.
!          IPNTR(2): POINTER TO THE CURRENT RESULT VECTOR Y.
!          IPNTR(3): POINTER TO THE VECTOR B * X WHEN USED IN THE
!                    SHIFT-AND-INVERT MODE.  X IS THE CURRENT OPERAND.
!          -------------------------------------------------------------
!
!  WORKD   COMPLEX*16 WORK ARRAY OF LENGTH 3*N.  (REVERSE COMMUNICATION)
!          DISTRIBUTED ARRAY TO BE USED IN THE BASIC ARNOLDI ITERATION
!          FOR REVERSE COMMUNICATION.  THE CALLING PROGRAM SHOULD NOT
!          USE WORKD AS TEMPORARY WORKSPACE DURING THE ITERATION !!!!!!
!          ON INPUT, WORKD(1:N) = B*RESID AND IS USED TO SAVE SOME
!          COMPUTATION AT THE FIRST STEP.
!
!  INFO    INTEGER.  (OUTPUT)
!          = 0: NORMAL EXIT.
!          > 0: SIZE OF THE SPANNING INVARIANT SUBSPACE OF OP FOUND.
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
!  2. R.B. LEHOUCQ, "ANALYSIS AND IMPLEMENTATION OF AN IMPLICITLY
!     RESTARTED ARNOLDI ITERATION", RICE UNIVERSITY TECHNICAL REPORT
!     TR95-13, DEPARTMENT OF COMPUTATIONAL AND APPLIED MATHEMATICS.
!
!\ROUTINES CALLED:
!     ZGETV0  ARPACK ROUTINE TO GENERATE THE INITIAL VECTOR.
!     IVOUT   ARPACK UTILITY ROUTINE THAT PRINTS INTEGERS.
!     ZMOUT   ARPACK UTILITY ROUTINE THAT PRINTS MATRICES
!     ZVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!       LAPACK ROUTINE THAT COMPUTES VARIOUS NORMS OF A MATRIX.
!     ZLASCL  LAPACK ROUTINE FOR CAREFUL SCALING OF A MATRIX.
!     DLABAD  LAPACK ROUTINE FOR DEFINING THE UNDERFLOW AND OVERFLOW
!             LIMITS.
!     DLAPY2  LAPACK ROUTINE TO COMPUTE SQRT(X**2+Y**2) CAREFULLY.
!     ZGEMV   LEVEL 2 BLAS ROUTINE FOR MATRIX VECTOR MULTIPLICATION.
!     ZAXPY   LEVEL 1 BLAS THAT COMPUTES A VECTOR TRIAD.
!     ZCOPY   LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER .
!     ZDOTC   LEVEL 1 BLAS THAT COMPUTES THE SCALAR PRODUCT OF TWO
!              VECTORS.
!     ZSCAL   LEVEL 1 BLAS THAT SCALES A VECTOR.
!     ZDSCAL  LEVEL 1 BLAS THAT SCALES A COMPLEX VECTOR BY A REAL
!             NUMBER.
!     DZNRM2  LEVEL 1 BLAS THAT COMPUTES THE NORM OF A VECTOR.
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
! FILE: NAITR.F   SID: 2.3   DATE OF SID: 8/27/96   RELEASE: 2
!
!\REMARKS
!  THE ALGORITHM IMPLEMENTED IS:
!
!  RESTART = .FALSE.
!  GIVEN V_{K} = [V_{1}, ..., V_{K}], R_{K};
!  R_{K} CONTAINS THE INITIAL RESIDUAL VECTOR EVEN FOR K = 0;
!  ALSO ASSUME THAT RNORM = || B*R_{K} || AND B*R_{K} ARE ALREADY
!  COMPUTED BY THE CALLING PROGRAM.
!
!  BETAJ = RNORM ; P_{K+1} = B*R_{K} ;
!  FOR  J = K+1, ..., K+NP  DO
!     1) IF ( BETAJ < TOL ) STOP OR RESTART DEPENDING ON J.
!        ( AT PRESENT TOL IS ZERO )
!        IF ( RESTART ) GENERATE A NEW STARTING VECTOR.
!     2) V_{J} = R(J-1)/BETAJ;  V_{J} = [V_{J-1}, V_{J}];
!        P_{J} = P_{J}/BETAJ
!     3) R_{J} = OP*V_{J} WHERE OP IS DEFINED AS IN ZNAUPD
!        FOR SHIFT-INVERT MODE P_{J} = B*V_{J} IS ALREADY AVAILABLE.
!        WNORM = || OP*V_{J} ||
!     4) COMPUTE THE J-TH STEP RESIDUAL VECTOR.
!        W_{J} =  V_{J}^T * B * OP * V_{J}
!        R_{J} =  OP*V_{J} - V_{J} * W_{J}
subroutine znaitr(ido, bmat, n, k, np, &
                  resid, rnorm, v, ldv, h, &
                  ldh, ipntr, workd, info, alpha)
!        H(:,J) = W_{J};
!        H(J,J-1) = RNORM
!        RNORM = || R_(J) ||
!        IF (RNORM > 0.717*WNORM) ACCEPT STEP AND GO BACK TO 1)
!     5) RE-ORTHOGONALIZATION STEP:
!        S = V_{J}'*B*R_{J}
!        R_{J} = R_{J} - V_{J}*S;  RNORM1 = || R_{J} ||
!        ALPHAJ = ALPHAJ + S_{J};
!     6) ITERATIVE REFINEMENT STEP:
!        IF (RNORM1 > 0.717*RNORM) THEN
!           RNORM = RNORM1
!           ACCEPT STEP AND GO BACK TO 1)
!        ELSE
!           RNORM = RNORM1
!           IF THIS IS THE FIRST TIME IN STEP 6), GOTO 5)
!           ELSE R_{J} LIES IN THE SPAN OF V_{J} NUMERICALLY.
!              SET R_{J} = 0 AND RNORM = 0; GOTO 1)
!        ENDIF
!  END DO
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!
!     %----------------------------------------------------%
!     | INCLUDE FILES FOR DEBUGGING AND TIMING INFORMATION |
!     %----------------------------------------------------%
!
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/dvout.h"
#include "asterfort/ivout.h"
#include "asterfort/zgetv0.h"
#include "asterfort/zmout.h"
#include "asterfort/zvout.h"
#include "blas/dlapy2.h"
#include "blas/dznrm2.h"
#include "blas/zaxpy.h"
#include "blas/zcopy.h"
#include "blas/zdotc.h"
#include "blas/zdscal.h"
#include "blas/zgemv.h"
#include "blas/zlanhs.h"
#include "blas/zlascl.h"
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
    integer(kind=8) :: ido, info, k, ldh, ldv, n, np
    real(kind=8) :: rnorm, alpha
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    integer(kind=8) :: ipntr(3)
    complex(kind=8) :: h(ldh, k+np), resid(n), v(ldv, k+np), workd(3*n)
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    complex(kind=8) :: one, zero
    real(kind=8) :: rone, rzero
    parameter(one=(1.0d+0, 0.0d+0), zero=(0.0d+0, 0.0d+0),&
     &           rone=1.0d+0, rzero=0.0d+0)
!
!     %--------------%
!     | LOCAL ARRAYS |
!     %--------------%
!
    real(kind=8) :: rtemp(2)
!
!     %---------------%
!     | LOCAL SCALARS |
!     %---------------%
!
    aster_logical :: first, orth1, orth2, rstart, step3, step4
    integer(kind=4) :: infol4
    integer(kind=8) :: ierr, i, ipj, irj, ivj, iter, itry, j, msglvl, jj
    real(kind=8) :: smlnum, tst1, ulp, unfl, betaj, temp1, rnorm1, wnorm, rbid(1)
    complex(kind=8) :: cnorm
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda, b_m
    blas_int :: b_kl, b_ku
!
    save first, orth1, orth2, rstart, step3, step4,&
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl,&
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
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
!        %-----------------------------------------%
!        | SET MACHINE-DEPENDENT CONSTANTS FOR THE |
!        | THE SPLITTING AND DEFLATION CRITERION.  |
!        | IF NORM(H) <= SQRT(OVFL),               |
!        | OVERFLOW SHOULD NOT OCCUR.              |
!        | REFERENCE: LAPACK SUBROUTINE ZLAHQR     |
!        %-----------------------------------------%
!
        unfl = r8miem()
! DUE TO CRS512         OVFL = ONE / UNFL
        ulp = r8prem()*0.5d0*isbaem()
        smlnum = unfl*(n/ulp)
        first = .false.
    end if
!
    if (ido .eq. 0) then
!
!        %-------------------------------%
!        | INITIALIZE TIMING STATISTICS  |
!        | & MESSAGE LEVEL FOR DEBUGGING |
!        %-------------------------------%
!
        msglvl = mnaitr
!
!        %------------------------------%
!        | INITIAL CALL TO THIS ROUTINE |
!        %------------------------------%
!
        info = 0
        step3 = .false.
        step4 = .false.
        rstart = .false.
        orth1 = .false.
        orth2 = .false.
        j = k+1
        ipj = 1
        irj = ipj+n
        ivj = irj+n
    end if
!
!     %-------------------------------------------------%
!     | WHEN IN REVERSE COMMUNICATION MODE ONE OF:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | WILL BE .TRUE. WHEN ....                        |
!     | STEP3: RETURN FROM COMPUTING OP*V_{J}.          |
!     | STEP4: RETURN FROM COMPUTING B-NORM OF OP*V_{J} |
!     | ORTH1: RETURN FROM COMPUTING B-NORM OF R_{J+1}  |
!     | ORTH2: RETURN FROM COMPUTING B-NORM OF          |
!     |        CORRECTION TO THE RESIDUAL VECTOR.       |
!     | RSTART: RETURN FROM OP COMPUTATIONS NEEDED BY   |
!     |         ZGETV0.                                 |
!     %-------------------------------------------------%
!
    if (step3) goto 50
    if (step4) goto 60
    if (orth1) goto 70
    if (orth2) goto 90
    if (rstart) goto 30
!
!     %-----------------------------%
!     | ELSE THIS IS THE FIRST STEP |
!     %-----------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | NOTE:  B*R_{J-1} IS ALREADY IN WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%
!
1000 continue
!
    if (msglvl .gt. 1) then
        call ivout(logfil, 1, [j], ndigit, '_NAITR: GENERATING ARNOLDI VECTOR NUMBER')
        call dvout(logfil, 1, [rnorm], ndigit, '_NAITR: B-NORM OF THE CURRENT RESIDUAL IS')
    end if
!
!        %---------------------------------------------------%
!        | STEP 1: CHECK IF THE B NORM OF J-TH RESIDUAL      |
!        | VECTOR IS ZERO. EQUIVALENT TO DETERMINE WHETHER   |
!        | AN EXACT J-STEP ARNOLDI FACTORIZATION IS PRESENT. |
!        %---------------------------------------------------%
!
    betaj = rnorm
    if (rnorm .gt. rzero) goto 40
!
!           %---------------------------------------------------%
!           | INVARIANT SUBSPACE FOUND, GENERATE A NEW STARTING |
!           | VECTOR WHICH IS ORTHOGONAL TO THE CURRENT ARNOLDI |
!           | BASIS AND CONTINUE THE ITERATION.                 |
!           %---------------------------------------------------%
!
    if (msglvl .gt. 0) then
        call ivout(logfil, 1, [j], ndigit, '_NAITR: ****** RESTART AT STEP ******')
    end if
!
!           %---------------------------------------------%
!           | ITRY IS THE LOOP VARIABLE THAT CONTROLS THE |
!           | MAXIMUM AMOUNT OF TIMES THAT A RESTART IS   |
!           | ATTEMPTED. NRSTRT IS USED BY STAT.H         |
!           %---------------------------------------------%
!
    betaj = rzero
    nrstrt = nrstrt+1
    itry = 1
20  continue
    rstart = .true.
    ido = 0
30  continue
!
!           %--------------------------------------%
!           | IF IN REVERSE COMMUNICATION MODE AND |
!           | RSTART = .TRUE. FLOW RETURNS HERE.   |
!           %--------------------------------------%
!
    call zgetv0(ido, bmat, .false._1, n, j, &
                v, ldv, resid, rnorm, ipntr, &
                workd, ierr, alpha)
    if (ido .ne. 99) goto 9000
    if (ierr .lt. 0) then
        itry = itry+1
        if (itry .le. 3) goto 20
!
!              %------------------------------------------------%
!              | GIVE UP AFTER SEVERAL RESTART ATTEMPTS.        |
!              | SET INFO TO THE SIZE OF THE INVARIANT SUBSPACE |
!              | WHICH SPANS OP AND EXIT.                       |
!              %------------------------------------------------%
!
        info = j-1
        ido = 99
        goto 9000
    end if
!
40  continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  V_{J} = R_{J-1}/RNORM AND P_{J} = P_{J}/RNORM  |
!        | NOTE THAT P_{J} = B*R_{J-1}. IN ORDER TO AVOID OVERFLOW |
!        | WHEN RECIPROCATING A SMALL RNORM, TEST AGAINST LOWER    |
!        | MACHINE BOUND.                                          |
!        %---------------------------------------------------------%
!
    b_n = to_blas_int(n)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, resid, b_incx, v(1, j), b_incy)
    if (rnorm .ge. unfl) then
        temp1 = rone/rnorm
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        call zdscal(b_n, temp1, v(1, j), b_incx)
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        call zdscal(b_n, temp1, workd(ipj), b_incx)
    else
!
!            %-----------------------------------------%
!            | TO SCALE BOTH V_{J} AND P_{J} CAREFULLY |
!            | USE LAPACK ROUTINE ZLASCL               |
!            %-----------------------------------------%
!
        b_lda = to_blas_int(n)
        b_kl = to_blas_int(1)
        b_ku = to_blas_int(1)
        b_m = to_blas_int(n)
        b_n = to_blas_int(1)
        call zlascl('G', b_kl, b_ku, rnorm, rone, &
                    b_m, b_n, v(1, j), b_lda, infol4)
        b_lda = to_blas_int(n)
        b_kl = to_blas_int(1)
        b_ku = to_blas_int(1)
        b_m = to_blas_int(n)
        b_n = to_blas_int(1)
        call zlascl('G', b_kl, b_ku, rnorm, rone, &
                    b_m, b_n, workd(ipj), b_lda, infol4)
    end if
!
!        %------------------------------------------------------%
!        | STEP 3:  R_{J} = OP*V_{J}; NOTE THAT P_{J} = B*V_{J} |
!        | NOTE THAT THIS IS NOT QUITE YET R_{J}. SEE STEP 4    |
!        %------------------------------------------------------%
!
    step3 = .true.
    nopx = nopx+1
    b_n = to_blas_int(n)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, v(1, j), b_incx, workd(ivj), b_incy)
    ipntr(1) = ivj
    ipntr(2) = irj
    ipntr(3) = ipj
    ido = 1
!
!        %-----------------------------------%
!        | EXIT IN ORDER TO COMPUTE OP*V_{J} |
!        %-----------------------------------%
!
    goto 9000
50  continue
!
!        %----------------------------------%
!        | BACK FROM REVERSE COMMUNICATION; |
!        | WORKD(IRJ:IRJ+N-1) := OP*V_{J}   |
!        | IF STEP3 = .TRUE.                |
!        %----------------------------------%
!
    step3 = .false.
!
!        %------------------------------------------%
!        | PUT ANOTHER COPY OF OP*V_{J} INTO RESID. |
!        %------------------------------------------%
!
    b_n = to_blas_int(n)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, workd(irj), b_incx, resid, b_incy)
!
!        %---------------------------------------%
!        | STEP 4:  FINISH EXTENDING THE ARNOLDI |
!        |          FACTORIZATION TO LENGTH J.   |
!        %---------------------------------------%
!
    if (bmat .eq. 'G') then
        nbx = nbx+1
        step4 = .true.
        ipntr(1) = irj
        ipntr(2) = ipj
        ido = 2
!
!           %-------------------------------------%
!           | EXIT IN ORDER TO COMPUTE B*OP*V_{J} |
!           %-------------------------------------%
!
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(ipj), b_incy)
    end if
60  continue
!
!        %----------------------------------%
!        | BACK FROM REVERSE COMMUNICATION; |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*V_{J} |
!        | IF STEP4 = .TRUE.                |
!        %----------------------------------%
!
    step4 = .false.
!
!        %-------------------------------------%
!        | THE FOLLOWING IS NEEDED FOR STEP 5. |
!        | COMPUTE THE B-NORM OF OP*V_{J}.     |
!        %-------------------------------------%
!
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cnorm = zdotc(b_n, resid, b_incx, workd(ipj), b_incy)
        wnorm = sqrt(dlapy2(dble(cnorm), dimag(cnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        wnorm = dznrm2(b_n, resid, b_incx)
    end if
!
!        %-----------------------------------------%
!        | COMPUTE THE J-TH RESIDUAL CORRESPONDING |
!        | TO THE J STEP FACTORIZATION.            |
!        | USE CLASSICAL GRAM SCHMIDT AND COMPUTE: |
!        | W_{J} <-  V_{J}^T * B * OP * V_{J}      |
!        | R_{J} <-  OP*V_{J} - V_{J} * W_{J}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | COMPUTE THE J FOURIER COEFFICIENTS W_{J} |
!        | WORKD(IPJ:IPJ+N-1) CONTAINS B*OP*V_{J}.  |
!        %------------------------------------------%
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('C', b_m, b_n, one, v, &
               b_lda, workd(ipj), b_incx, zero, h(1, j), &
               b_incy)
!
!        %--------------------------------------%
!        | ORTHOGONALIZE R_{J} AGAINST V_{J}.   |
!        | RESID CONTAINS OP*V_{J}. SEE STEP 3. |
!        %--------------------------------------%
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('N', b_m, b_n, -one, v, &
               b_lda, h(1, j), b_incx, one, resid, &
               b_incy)
!
    if (j .gt. 1) h(j, j-1) = dcmplx(betaj, rzero)
!
!
    orth1 = .true.
!
    if (bmat .eq. 'G') then
        nbx = nbx+1
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(irj), b_incy)
        ipntr(1) = irj
        ipntr(2) = ipj
        ido = 2
!
!           %----------------------------------%
!           | EXIT IN ORDER TO COMPUTE B*R_{J} |
!           %----------------------------------%
!
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(ipj), b_incy)
    end if
70  continue
!
!        %---------------------------------------------------%
!        | BACK FROM REVERSE COMMUNICATION IF ORTH1 = .TRUE. |
!        | WORKD(IPJ:IPJ+N-1) := B*R_{J}.                    |
!        %---------------------------------------------------%
!
    orth1 = .false.
!
!        %------------------------------%
!        | COMPUTE THE B-NORM OF R_{J}. |
!        %------------------------------%
!
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cnorm = zdotc(b_n, resid, b_incx, workd(ipj), b_incy)
        rnorm = sqrt(dlapy2(dble(cnorm), dimag(cnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        rnorm = dznrm2(b_n, resid, b_incx)
    end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: RE-ORTHOGONALIZATION / ITERATIVE REFINEMENT PHASE |
!        | MAXIMUM NITER_ITREF TRIES.                                |
!        |                                                           |
!        |          S      = V_{J}^T * B * R_{J}                     |
!        |          R_{J}  = R_{J} - V_{J}*S                         |
!        |          ALPHAJ = ALPHAJ + S_{J}                          |
!        |                                                           |
!        | THE STOPPING CRITERIA USED FOR ITERATIVE REFINEMENT IS    |
!        | DISCUSSED IN PARLETT'S BOOK SEP, PAGE 107 AND IN GRAGG &  |
!        | REICHEL ACM TOMS PAPER; ALGORITHM 686, DEC. 1990.         |
!        | DETERMINE IF WE NEED TO CORRECT THE RESIDUAL. THE GOAL IS |
!        | TO ENFORCE ||V(:,1:J)^T * R_{J}|| .LE. EPS * || R_{J} ||  |
!        | THE FOLLOWING TEST DETERMINES WHETHER THE SINE OF THE     |
!        | ANGLE BETWEEN  OP*X AND THE COMPUTED RESIDUAL IS LESS     |
!        | THAN OR EQUAL TO 0.717.                                   |
!        %-----------------------------------------------------------%
!
    if (rnorm .gt. alpha*wnorm) goto 100
!
    iter = 0
    nrorth = nrorth+1
!
!        %---------------------------------------------------%
!        | ENTER THE ITERATIVE REFINEMENT PHASE. IF FURTHER  |
!        | REFINEMENT IS NECESSARY, LOOP BACK HERE. THE LOOP |
!        | VARIABLE IS ITER. PERFORM A STEP OF CLASSICAL     |
!        | GRAM-SCHMIDT USING ALL THE ARNOLDI VECTORS V_{J}  |
!        %---------------------------------------------------%
!
80  continue
!
    if (msglvl .gt. 2) then
        rtemp(1) = wnorm
        rtemp(2) = rnorm
        call dvout(logfil, 2, rtemp, ndigit, '_NAITR: RE-ORTHOGONALIZATION; WNORM AND RNORM ARE')
        call zvout(logfil, j, h(1, j), ndigit, '_NAITR: J-TH COLUMN OF H')
    end if
!
!        %----------------------------------------------------%
!        | COMPUTE V_{J}^T * B * R_{J}.                       |
!        | WORKD(IRJ:IRJ+J-1) = V(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('C', b_m, b_n, one, v, &
               b_lda, workd(ipj), b_incx, zero, workd(irj), &
               b_incy)
!
!        %---------------------------------------------%
!        | COMPUTE THE CORRECTION TO THE RESIDUAL:     |
!        | R_{J} = R_{J} - V_{J} * WORKD(IRJ:IRJ+J-1). |
!        | THE CORRECTION TO H IS V(:,1:J)*H(1:J,1:J)  |
!        | + V(:,1:J)*WORKD(IRJ:IRJ+J-1)*E'_J.         |
!        %---------------------------------------------%
!
    b_lda = to_blas_int(ldv)
    b_m = to_blas_int(n)
    b_n = to_blas_int(j)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zgemv('N', b_m, b_n, -one, v, &
               b_lda, workd(irj), b_incx, one, resid, &
               b_incy)
    b_n = to_blas_int(j)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zaxpy(b_n, one, workd(irj), b_incx, h(1, j), &
               b_incy)
!
    orth2 = .true.
    if (bmat .eq. 'G') then
        nbx = nbx+1
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(irj), b_incy)
        ipntr(1) = irj
        ipntr(2) = ipj
        ido = 2
!
!           %-----------------------------------%
!           | EXIT IN ORDER TO COMPUTE B*R_{J}. |
!           | R_{J} IS THE CORRECTED RESIDUAL.  |
!           %-----------------------------------%
!
        goto 9000
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, resid, b_incx, workd(ipj), b_incy)
    end if
90  continue
!
!        %---------------------------------------------------%
!        | BACK FROM REVERSE COMMUNICATION IF ORTH2 = .TRUE. |
!        %---------------------------------------------------%
!
!        %-----------------------------------------------------%
!        | COMPUTE THE B-NORM OF THE CORRECTED RESIDUAL R_{J}. |
!        %-----------------------------------------------------%
!
    if (bmat .eq. 'G') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cnorm = zdotc(b_n, resid, b_incx, workd(ipj), b_incy)
        rnorm1 = sqrt(dlapy2(dble(cnorm), dimag(cnorm)))
    else if (bmat .eq. 'I') then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        rnorm1 = dznrm2(b_n, resid, b_incx)
    end if
!
    if (msglvl .gt. 0 .and. iter .gt. 0) then
        call ivout(logfil, 1, [j], ndigit, '_NAITR: ITERATIVE REFINEMENT FOR ARNOLDI RESIDUAL')
        if (msglvl .gt. 2) then
            rtemp(1) = rnorm
            rtemp(2) = rnorm1
            call dvout(logfil, 2, rtemp, ndigit, &
                       '_NAITR: ITERATIVE REFINEMENT ; RNORM AND RNORM1 ARE')
        end if
    end if
!
!        %-----------------------------------------%
!        | DETERMINE IF WE NEED TO PERFORM ANOTHER |
!        | STEP OF RE-ORTHOGONALIZATION.           |
!        %-----------------------------------------%
!
    if (rnorm1 .gt. alpha*rnorm) then
!
!           %---------------------------------------%
!           | NO NEED FOR FURTHER REFINEMENT.       |
!           | THE COSINE OF THE ANGLE BETWEEN THE   |
!           | CORRECTED RESIDUAL VECTOR AND THE OLD |
!           | RESIDUAL VECTOR IS GREATER THAN 0.717 |
!           | IN OTHER WORDS THE CORRECTED RESIDUAL |
!           | AND THE OLD RESIDUAL VECTOR SHARE AN  |
!           | ANGLE OF LESS THAN ARCCOS(0.717)      |
!           %---------------------------------------%
!
        rnorm = rnorm1
!
    else
!
!           %-------------------------------------------%
!           | ANOTHER STEP OF ITERATIVE REFINEMENT STEP |
!           | IS REQUIRED. NITREF IS USED BY STAT.H     |
!           %-------------------------------------------%
!
        nitref = nitref+1
        rnorm = rnorm1
        iter = iter+1
        if (iter .le. 1) goto 80
!
!           %-------------------------------------------------%
!           | OTHERWISE RESID IS NUMERICALLY IN THE SPAN OF V |
!           %-------------------------------------------------%
!
        do jj = 1, n
            resid(jj) = zero
        end do
        rnorm = rzero
    end if
!
!        %----------------------------------------------%
!        | BRANCH HERE DIRECTLY IF ITERATIVE REFINEMENT |
!        | WASN'T NECESSARY OR AFTER AT MOST NITER_REF  |
!        | STEPS OF ITERATIVE REFINEMENT.               |
!        %----------------------------------------------%
!
100 continue
!
    rstart = .false.
    orth2 = .false.
!
!        %------------------------------------%
!        | STEP 6: UPDATE  J = J+1;  CONTINUE |
!        %------------------------------------%
!
    j = j+1
    if (j .gt. k+np) then
        ido = 99
        do i = max(1, k), k+np-1
!
!              %--------------------------------------------%
!              | CHECK FOR SPLITTING AND DEFLATION.         |
!              | USE A STANDARD TEST AS IN THE QR ALGORITHM |
!              | REFERENCE: LAPACK SUBROUTINE ZLAHQR        |
!              %--------------------------------------------%
!
            tst1 = dlapy2( &
                   dble(h(i, i)), dimag(h(i, i)))+dlapy2(dble(h(i+1, i+1)), dimag(h(i+1, i+1)))
            b_lda = to_blas_int(ldh)
            b_n = to_blas_int(k+np)
            if (tst1 .eq. dble(zero)) tst1 = zlanhs('1', b_n, h, b_lda, rbid)
            if (dlapy2(dble(h(i+1, i)), dimag(h(i+1, i))) .le. max(ulp*tst1, smlnum)) h(i+1, i) = &
                zero
        end do
!
        if (msglvl .gt. 2) then
            call zmout(logfil, k+np, k+np, h, ldh, &
                       ndigit, '_NAITR: FINAL UPPER HESSENBERG MATRIX H OF ORDER K+NP')
        end if
!
        goto 9000
    end if
!
!        %--------------------------------------------------------%
!        | LOOP BACK TO EXTEND THE FACTORIZATION BY ANOTHER STEP. |
!        %--------------------------------------------------------%
!
    goto 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
9000 continue
    call matfpe(1)
!
!     %---------------%
!     | END OF ZNAITR |
!     %---------------%
!
end subroutine
