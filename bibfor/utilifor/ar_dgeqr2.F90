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
subroutine ar_dgeqr2(m, n, a, lda, tau, &
                     work, info)
!
!     SUBROUTINE LAPACK CALCULANT UNE FACTORISATION QR.
!-----------------------------------------------------------------------
!  -- LAPACK ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992
!
!  PURPOSE
!  =======
!
!  DGEQR2 COMPUTES A QR FACTORIZATION OF A REAL M BY N MATRIX A:
!  A = Q * R.
!
!  ARGUMENTS
!  =========
!
!  M       (INPUT) INTEGER
!          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
!
!  N       (INPUT) INTEGER
!          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
!
!  A       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDA,N)
!          ON ENTRY, THE M BY N MATRIX A.
!          ON EXIT, THE ELEMENTS ON AND ABOVE THE DIAGONAL OF THE ARRAY
!          CONTAIN THE MIN(M,N) BY N UPPER TRAPEZOIDAL MATRIX R (R IS
!          UPPER TRIANGULAR IF M >= N), THE ELEMENTS BELOW THE DIAGONAL,
!          WITH THE ARRAY TAU, REPRESENT THE ORTHOGONAL MATRIX Q AS A
!          PRODUCT OF ELEMENTARY REFLECTORS (SEE FURTHER DETAILS).
!
!  LDA     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
!
!  TAU     (OUTPUT) REAL*8 ARRAY, DIMENSION (MIN(M,N))
!          THE SCALAR FACTORS OF THE ELEMENTARY REFLECTORS (SEE FURTHER
!          DETAILS).
!
!  WORK    (WORKSPACE) REAL*8 ARRAY, DIMENSION (N)
!
!  INFO    (OUTPUT) INTEGER
!          = 0: SUCCESSFUL EXIT
!          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
!
!  FURTHER DETAILS
!  ===============
!
!  THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY REFLECTORS
!
!     Q = H(1) H(2) . . . H(K), WHERE K = MIN(M,N).
!
!  EACH H(I) HAS THE FORM
!
!     H(I) = I - TAU * V * V'
!
!  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
!  V(1:I-1) = 0 AND V(I) = 1, V(I+1:M) IS STORED ON EXIT IN A(I+1:M,I),
!  AND TAU IN TAU(I).
!
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 1 RETURN PAR GOTO 1000,
!            IMPLICIT NONE.
! INTRINSIC FUNCTIONS
!            MAX, MIN.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterc/matfpe.h"
#include "asterfort/ar_dlarfg.h"
#include "asterfort/xerbla.h"
#include "blas/dlarf.h"
    integer(kind=8) :: info, lda, m, n
!     ..
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: a(lda, *), tau(*), work(*)
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: one
    parameter(one=1.0d+0)
!     ..
!     .. LOCAL SCALARS ..
    integer(kind=8) :: i, k
    real(kind=8) :: aii
    blas_int :: b_incv, b_ldc, b_m, b_n
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
!     TEST THE INPUT ARGUMENTS
!
    info = 0
    if (m .lt. 0) then
        info = -1
    else if (n .lt. 0) then
        info = -2
    else if (lda .lt. max(1, m)) then
        info = -4
    end if
    if (info .ne. 0) then
        call xerbla('DGEQR2', -info)
        goto 1000
    end if
!
    k = min(m, n)
!
    do i = 1, k
!
!        GENERATE ELEMENTARY REFLECTOR H(I) TO ANNIHILATE A(I+1:M,I)
!
        call ar_dlarfg(m-i+1, a(i, i), a(min(i+1, m), i), 1, tau(i))
        if (i .lt. n) then
!
!           APPLY H(I) TO A(I:M,I+1:N) FROM THE LEFT
!
            aii = a(i, i)
            a(i, i) = one
            b_ldc = to_blas_int(lda)
            b_m = to_blas_int(m-i+1)
            b_n = to_blas_int(n-i)
            b_incv = to_blas_int(1)
            call dlarf('L', b_m, b_n, a(i, i), b_incv, &
                       tau(i), a(i, i+1), b_ldc, work)
            a(i, i) = aii
        end if
    end do
1000 continue
!
    call matfpe(1)
!
!     END OF DGEQR2
!
end subroutine
