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
subroutine ar_ztrexc(compq, n, t, ldt, q, &
                     ldq, ifst, ilst, info)
!  -- LAPACK ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     MARCH 31, 1993
!
!  PURPOSE
!  =======
!
!  ZTREXC REORDERS THE SCHUR FACTORIZATION OF A COMPLEX MATRIX
!  A = Q*T*Q**H, SO THAT THE DIAGONAL ELEMENT OF T WITH ROW INDEX IFST
!  IS MOVED TO ROW ILST.
!
!  THE SCHUR FORM T IS REORDERED BY A UNITARY SIMILARITY TRANSFORMATION
!  Z**H*T*Z, AND OPTIONALLY THE MATRIX Q OF SCHUR VECTORS IS UPDATED BY
!  POSTMULTPLYING IT WITH Z.
!
!  ARGUMENTS
!  =========
!
!  COMPQ   (INPUT) CHARACTER*1
!          = 'V':  UPDATE THE MATRIX Q OF SCHUR VECTORS;
!          = 'N':  DO NOT UPDATE Q.
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX T. N >= 0.
!
!  T       (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDT,N)
!          ON ENTRY, THE UPPER TRIANGULAR MATRIX T.
!          ON EXIT, THE REORDERED UPPER TRIANGULAR MATRIX.
!
!  LDT     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY T. LDT >= MAX(1,N).
!
!  Q       (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDQ,N)
!          ON ENTRY, IF COMPQ = 'V', THE MATRIX Q OF SCHUR VECTORS.
!          ON EXIT, IF COMPQ = 'V', Q HAS BEEN POSTMULTIPLIED BY THE
!          UNITARY TRANSFORMATION MATRIX Z WHICH REORDERS T.
!          IF COMPQ = 'N', Q IS NOT REFERENCED.
!
!  LDQ     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY Q.  LDQ >= MAX(1,N).
!
!  IFST    (INPUT) INTEGER
!  ILST    (INPUT) INTEGER
!          SPECIFY THE REORDERING OF THE DIAGONAL ELEMENTS OF T:
!          THE ELEMENT WITH ROW INDEX IFST IS MOVED TO ROW ILST BY A
!          SEQUENCE OF TRANSPOSITIONS BETWEEN ADJACENT ELEMENTS.
!          1 <= IFST <= N; 1 <= ILST <= N.
!
!  INFO    (OUTPUT) INTEGER
!          = 0:  SUCCESSFUL EXIT
!          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
!
!  =====================================================================
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 1 RETURN PAR GOTO 1000,
!            IMPLICIT NONE.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/ar_zlartg.h"
#include "asterfort/xerbla.h"
#include "blas/lsame.h"
#include "blas/zrot.h"
    character(len=1) :: compq
    integer(kind=8) :: ifst, ilst, info, ldq, ldt, n
!     ..
!     .. ARRAY ARGUMENTS ..
    complex(kind=8) :: q(ldq, *), t(ldt, *)
!     ..
!
!     .. LOCAL SCALARS ..
    aster_logical :: wantq
    integer(kind=8) :: k, m1, m2, m3
    real(kind=8) :: cs
    complex(kind=8) :: sn, t11, t22, temp
    blas_int :: b_incx, b_incy, b_n
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
!     DECODE AND TEST THE INPUT PARAMETERS.
!
    info = 0
    wantq = lsame(compq, 'V')
    if (.not. lsame(compq, 'N') .and. .not. wantq) then
        info = -1
    else if (n .lt. 0) then
        info = -2
    else if (ldt .lt. max(1, n)) then
        info = -4
    else if (ldq .lt. 1 .or. (wantq .and. ldq .lt. max(1, n))) then
        info = -6
    else if (ifst .lt. 1 .or. ifst .gt. n) then
        info = -7
    else if (ilst .lt. 1 .or. ilst .gt. n) then
        info = -8
    end if
    if (info .ne. 0) then
        call xerbla('ZTREXC', -info)
        goto 1000
    end if
!
!     QUICK RETURN IF POSSIBLE
!
    if (n .eq. 1 .or. ifst .eq. ilst) goto 1000
!
    if (ifst .lt. ilst) then
!
!        MOVE THE IFST-TH DIAGONAL ELEMENT FORWARD DOWN THE DIAGONAL.
!
        m1 = 0
        m2 = -1
        m3 = 1
    else
!
!        MOVE THE IFST-TH DIAGONAL ELEMENT BACKWARD UP THE DIAGONAL.
!
        m1 = -1
        m2 = 0
        m3 = -1
    end if
!
    do k = ifst+m1, ilst+m2, m3
!
!        INTERCHANGE THE K-TH AND (K+1)-TH DIAGONAL ELEMENTS.
!
        t11 = t(k, k)
        t22 = t(k+1, k+1)
!
!        DETERMINE THE TRANSFORMATION TO PERFORM THE INTERCHANGE.
!
        call ar_zlartg(t(k, k+1), t22-t11, cs, sn, temp)
!
!        APPLY TRANSFORMATION TO THE MATRIX T.
!
        if (k+2 .le. n) then
            b_n = to_blas_int(n-k-1)
            b_incx = to_blas_int(ldt)
            b_incy = to_blas_int(ldt)
            call zrot(b_n, t(k, k+2), b_incx, t(k+1, k+2), b_incy, &
                      cs, sn)
        end if
        b_n = to_blas_int(k-1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zrot(b_n, t(1, k), b_incx, t(1, k+1), b_incy, &
                  cs, dconjg(sn))
!
        t(k, k) = t22
        t(k+1, k+1) = t11
!
        if (wantq) then
!
!           ACCUMULATE TRANSFORMATION IN THE MATRIX Q.
!
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zrot(b_n, q(1, k), b_incx, q(1, k+1), b_incy, &
                      cs, dconjg(sn))
        end if
!
    end do
!
1000 continue
    call matfpe(1)
!
!     END OF ZTREXC
!
end subroutine
