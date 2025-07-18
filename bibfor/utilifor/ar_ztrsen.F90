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
subroutine ar_ztrsen(select, n, t, ldt, q, &
                     ldq, w, m, info)
!  -- LAPACK ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     MARCH 31, 1993
!
!
!  PURPOSE
!  =======
!
!  ZTRSEN REORDERS THE SCHUR FACTORIZATION OF A COMPLEX MATRIX
!  A = Q*T*Q**H, SO THAT A SELECTED CLUSTER OF EIGENVALUES APPEARS IN
!  THE LEADING POSITIONS ON THE DIAGONAL OF THE UPPER TRIANGULAR MATRIX
!  T, AND THE LEADING COLUMNS OF Q FORM AN ORTHONORMAL BASIS OF THE
!  CORRESPONDING RIGHT INVARIANT SUBSPACE.
!
!  OPTIONALLY THE ROUTINE COMPUTES THE RECIPROCAL CONDITION NUMBERS OF
!  THE CLUSTER OF EIGENVALUES AND/OR THE INVARIANT SUBSPACE.
!
!  ARGUMENTS
!  =========
!
!  SELECT  (INPUT) LOGICAL ARRAY, DIMENSION (N)
!          SELECT SPECIFIES THE EIGENVALUES IN THE SELECTED CLUSTER. TO
!          SELECT THE J-TH EIGENVALUE, SELECT(J) MUST BE SET TO .TRUE..
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX T. N >= 0.
!
!  T       (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDT,N)
!          ON ENTRY, THE UPPER TRIANGULAR MATRIX T.
!          ON EXIT, T IS OVERWRITTEN BY THE REORDERED MATRIX T, WITH THE
!          SELECTED EIGENVALUES AS THE LEADING DIAGONAL ELEMENTS.
!
!  LDT     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY T. LDT >= MAX(1,N).
!
!  Q       (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDQ,N)
!          ON ENTRY, THE MATRIX Q OF SCHUR VECTORS.
!          ON EXIT, Q HAS BEEN POSTMULTIPLIED BY THE
!          UNITARY TRANSFORMATION MATRIX WHICH REORDERS T; THE LEADING M
!          COLUMNS OF Q FORM AN ORTHONORMAL BASIS FOR THE SPECIFIED
!          INVARIANT SUBSPACE.
!
!  LDQ     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY Q.
!          LDQ >= N.
!
!  W       (OUTPUT) COMPLEX*16
!          THE REORDERED EIGENVALUES OF T, IN THE SAME ORDER AS THEY
!          APPEAR ON THE DIAGONAL OF T.
!
!  M       (OUTPUT) INTEGER
!          THE DIMENSION OF THE SPECIFIED INVARIANT SUBSPACE.
!          0 <= M <= N.
!
!  INFO    (OUTPUT) INTEGER
!          = 0:  SUCCESSFUL EXIT
!          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
!
!  FURTHER DETAILS
!  ===============
!
!  ZTRSEN FIRST COLLECTS THE SELECTED EIGENVALUES BY COMPUTING A UNITARY
!  TRANSFORMATION Z TO MOVE THEM TO THE TOP LEFT CORNER OF T. IN OTHER
!  WORDS, THE SELECTED EIGENVALUES ARE THE EIGENVALUES OF T11 IN:
!
!                Z'*T*Z = ( T11 T12 ) N1
!                         (  0  T22 ) N2
!                            N1  N2
!
!  WHERE N = N1+N2 AND Z' MEANS THE CONJUGATE TRANSPOSE OF Z. THE FIRST
!  N1 COLUMNS OF Z SPAN THE SPECIFIED INVARIANT SUBSPACE OF T.
!
!  IF T HAS BEEN OBTAINED FROM THE SCHUR FACTORIZATION OF A MATRIX
!  A = Q*T*Q', THEN THE REORDERED SCHUR FACTORIZATION OF A IS GIVEN BY
!  A = (Q*Z)*(Z'*T*Z)*(Q*Z)', AND THE FIRST N1 COLUMNS OF Q*Z SPAN THE
!  CORRESPONDING INVARIANT SUBSPACE OF A.
!
!  THE RECIPROCAL CONDITION NUMBER OF THE AVERAGE OF THE EIGENVALUES OF
!  T11 MAY BE RETURNED IN S. S LIES BETWEEN 0 (VERY BADLY CONDITIONED)
!  AND 1 (VERY WELL CONDITIONED). IT IS COMPUTED AS FOLLOWS. FIRST WE
!  COMPUTE R SO THAT
!
!                         P = ( I  R ) N1
!                             ( 0  0 ) N2
!                               N1 N2
!
!  IS THE PROJECTOR ON THE INVARIANT SUBSPACE ASSOCIATED WITH T11.
!  R IS THE SOLUTION OF THE SYLVESTER EQUATION:
!
!                        T11*R - R*T22 = T12.
!
!  LET F-NORM(M) DENOTE THE FROBENIUS-NORM OF M AND 2-NORM(M) DENOTE
!  THE TWO-NORM OF M. THEN S IS COMPUTED AS THE LOWER BOUND
!
!                      (1 + F-NORM(R)**2)**(-1/2)
!
!  ON THE RECIPROCAL OF 2-NORM(P), THE TRUE RECIPROCAL CONDITION NUMBER.
!  S CANNOT UNDERESTIMATE 1 / 2-NORM(P) BY MORE THAN A FACTOR OF
!  SQRT(N).
!
!  AN APPROXIMATE ERROR BOUND FOR THE COMPUTED AVERAGE OF THE
!  EIGENVALUES OF T11 IS
!
!                         EPS * NORM(T) / S
!
!  WHERE EPS IS THE MACHINE PRECISION.
!
!  THE RECIPROCAL CONDITION NUMBER OF THE RIGHT INVARIANT SUBSPACE
!  SPANNED BY THE FIRST N1 COLUMNS OF Z (OR OF Q*Z) IS RETURNED IN SEP.
!  SEP IS DEFINED AS THE SEPARATION OF T11 AND T22:
!
!                     SEP( T11, T22 ) = SIGMA-MIN( C )
!
!  WHERE SIGMA-MIN(C) IS THE SMALLEST SINGULAR VALUE OF THE
!  N1*N2-BY-N1*N2 MATRIX
!
!     C  = KPROD( I(N2), T11 ) - KPROD( TRANSPOSE(T22), I(N1) )
!
!  I(M) IS AN M BY M IDENTITY MATRIX, AND KPROD DENOTES THE KRONECKER
!  PRODUCT. WE ESTIMATE SIGMA-MIN(C) BY THE RECIPROCAL OF AN ESTIMATE OF
!  THE 1-NORM OF INVERSE(C). THE TRUE RECIPROCAL 1-NORM OF INVERSE(C)
!  CANNOT DIFFER FROM SIGMA-MIN(C) BY MORE THAN A FACTOR OF SQRT(N1*N2).
!
!  WHEN SEP IS SMALL, SMALL CHANGES IN T CAN CAUSE LARGE CHANGES IN
!  THE INVARIANT SUBSPACE. AN APPROXIMATE BOUND ON THE MAXIMUM ANGULAR
!  ERROR IN THE COMPUTED RIGHT INVARIANT SUBSPACE IS
!
!                      EPS * NORM(T) / SEP
!
!  =====================================================================
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 1 RETURN PAR GOTO 100,
!            IMPLICIT NONE.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/ar_ztrexc.h"
#include "asterfort/xerbla.h"
    integer(kind=8) :: info, ldq, ldt, m, n
!     ..
!     .. ARRAY ARGUMENTS ..
    aster_logical :: select(*)
    complex(kind=8) :: q(ldq, *), t(ldt, *), w(*)
!     ..
!     .. LOCAL SCALARS ..
    integer(kind=8) :: ierr, k, ks, n1, n2, nn
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
!     DECODE AND TEST THE INPUT PARAMETERS.
!
!     SET M TO THE NUMBER OF SELECTED EIGENVALUES.
!
    m = 0
    do k = 1, n
        if (select(k)) m = m+1
    end do
!
    n1 = m
    n2 = n-m
    nn = n1*n2
!
    info = 0
    if (n .lt. 0) then
        info = -4
    else if (ldt .lt. max(1, n)) then
        info = -6
    else if (ldq .lt. 1 .or. (ldq .lt. n)) then
        info = -8
    end if
    if (info .ne. 0) then
        call xerbla('ZTRSEN', -info)
        goto 100
    end if
!
!     QUICK RETURN IF POSSIBLE
!
    if (m .eq. n .or. m .eq. 0) then
        goto 40
    end if
!
!     COLLECT THE SELECTED EIGENVALUES AT THE TOP LEFT CORNER OF T.
!
    ks = 0
    do k = 1, n
        if (select(k)) then
            ks = ks+1
!
!           SWAP THE K-TH EIGENVALUE TO POSITION KS.
!
            if (k .ne. ks) call ar_ztrexc('V', n, t, ldt, q, &
                                          ldq, k, ks, ierr)
        end if
    end do
!
!
40  continue
!
!     COPY REORDERED EIGENVALUES TO W.
!
    do k = 1, n
        w(k) = t(k, k)
    end do
100 continue
    call matfpe(1)
!
!     END OF ZTRSEN
!
end subroutine
