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
subroutine zsortc(which, apply, n, x, y)
!
!     SUBROUTINE ARPACK TRIANT DES VECTEURS COMPLEXES.
!-----------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZSORTC
!
!\DESCRIPTION:
!  SORTS THE COMPLEX*16 ARRAY IN X INTO THE ORDER
!  SPECIFIED BY WHICH AND OPTIONALLY APPLIES THE PERMUTATION TO THE
!  DOUBLE PRECISION  ARRAY Y.
!
!\USAGE:
!  CALL ZSORTC
!     ( WHICH, APPLY, N, X, Y )
!
!\ARGUMENTS
!  WHICH   CHARACTER*2.  (INPUT)
!          'LM' -> SORT X INTO INCREASING ORDER OF MAGNITUDE.
!          'SM' -> SORT X INTO DECREASING ORDER OF MAGNITUDE.
!          'LR' -> SORT X WITH REAL(X) IN INCREASING ALGEBRAIC ORDER
!          'SR' -> SORT X WITH REAL(X) IN DECREASING ALGEBRAIC ORDER
!          'LI' -> SORT X WITH IMAG(X) IN INCREASING ALGEBRAIC ORDER
!          'SI' -> SORT X WITH IMAG(X) IN DECREASING ALGEBRAIC ORDER
!
!  APPLY   LOGICAL.  (INPUT)
!          APPLY = .TRUE.  -> APPLY THE SORTED ORDER TO ARRAY Y.
!          APPLY = .FALSE. -> DO NOT APPLY THE SORTED ORDER TO ARRAY Y.
!
!  N       INTEGER.  (INPUT)
!          SIZE OF THE ARRAYS.
!
!  X       COMPLEX*16 ARRAY OF LENGTH N.  (INPUT/OUTPUT)
!          THIS IS THE ARRAY TO BE SORTED.
!
!  Y       COMPLEX*16 ARRAY OF LENGTH N.  (INPUT/OUTPUT)
!
!\ENDDOC
!
!-----------------------------------------------------------------------
!
!\BEGINLIB
!
!\ROUTINES CALLED:
!     DLAPY2  LAPACK ROUTINE TO COMPUTE SQRT(X**2+Y**2) CAREFULLY.
!
!\AUTHOR
!     DANNY SORENSEN               PHUONG VU
!     RICHARD LEHOUCQ              CRPC / RICE UNIVERSITY
!     DEPT. OF COMPUTATIONAL &     HOUSTON, TEXAS
!     APPLIED MATHEMATICS
!     RICE UNIVERSITY
!     HOUSTON, TEXAS
!
!     ADAPTED FROM THE SORT ROUTINE IN LANSO.
!
!\SCCS INFORMATION: @(#)
! FILE: SORTC.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "blas/dlapy2.h"
    character(len=2) :: which
    aster_logical :: apply
    integer(kind=8) :: n
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    complex(kind=8) :: x(0:n-1), y(0:n-1)
!
!     %---------------%
!     | LOCAL SCALARS |
!     %---------------%
!
    integer(kind=8) :: i, igap, j
    complex(kind=8) :: temp
    real(kind=8) :: temp1, temp2
!
!     %--------------------%
!     | EXTERNAL FUNCTIONS |
!     %--------------------%
!
!
!
!     %-----------------------%
!     | EXECUTABLE STATEMENTS |
!     %-----------------------%
!
    call matfpe(-1)
    igap = n/2
!
    if (which .eq. 'LM') then
!
!        %--------------------------------------------%
!        | SORT X INTO INCREASING ORDER OF MAGNITUDE. |
!        %--------------------------------------------%
!
10      continue
        if (igap .eq. 0) goto 9000
!
        do i = igap, n-1
            j = i-igap
20          continue
!
            if (j .lt. 0) goto 30
!
            temp1 = dlapy2(dble(x(j)), dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)), dimag(x(j+igap)))
!
            if (temp1 .gt. temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 30
            end if
            j = j-igap
            goto 20
30          continue
        end do
        igap = igap/2
        goto 10
!
    else if (which .eq. 'SM') then
!
!        %--------------------------------------------%
!        | SORT X INTO DECREASING ORDER OF MAGNITUDE. |
!        %--------------------------------------------%
!
40      continue
        if (igap .eq. 0) goto 9000
!
        do i = igap, n-1
            j = i-igap
50          continue
!
            if (j .lt. 0) goto 60
!
            temp1 = dlapy2(dble(x(j)), dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)), dimag(x(j+igap)))
!
            if (temp1 .lt. temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 60
            end if
            j = j-igap
            goto 50
60          continue
        end do
        igap = igap/2
        goto 40
!
    else if (which .eq. 'LR') then
!
!        %------------------------------------------------%
!        | SORT XREAL INTO INCREASING ORDER OF ALGEBRAIC. |
!        %------------------------------------------------%
!
70      continue
        if (igap .eq. 0) goto 9000
!
        do i = igap, n-1
            j = i-igap
80          continue
!
            if (j .lt. 0) goto 90
!
            if (dble(x(j)) .gt. dble(x(j+igap))) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 90
            end if
            j = j-igap
            goto 80
90          continue
        end do
        igap = igap/2
        goto 70
!
    else if (which .eq. 'SR') then
!
!        %------------------------------------------------%
!        | SORT XREAL INTO DECREASING ORDER OF ALGEBRAIC. |
!        %------------------------------------------------%
!
100     continue
        if (igap .eq. 0) goto 9000
        do i = igap, n-1
            j = i-igap
110         continue
!
            if (j .lt. 0) goto 120
!
            if (dble(x(j)) .lt. dble(x(j+igap))) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 120
            end if
            j = j-igap
            goto 110
120         continue
        end do
        igap = igap/2
        goto 100
!
    else if (which .eq. 'LI') then
!
!        %--------------------------------------------%
!        | SORT XIMAG INTO INCREASING ALGEBRAIC ORDER |
!        %--------------------------------------------%
!
130     continue
        if (igap .eq. 0) goto 9000
        do i = igap, n-1
            j = i-igap
140         continue
!
            if (j .lt. 0) goto 150
!
            if (dimag(x(j)) .gt. dimag(x(j+igap))) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 150
            end if
            j = j-igap
            goto 140
150         continue
        end do
        igap = igap/2
        goto 130
!
    else if (which .eq. 'SI') then
!
!        %---------------------------------------------%
!        | SORT XIMAG INTO DECREASING ALGEBRAIC ORDER  |
!        %---------------------------------------------%
!
160     continue
        if (igap .eq. 0) goto 9000
        do i = igap, n-1
            j = i-igap
170         continue
!
            if (j .lt. 0) goto 180
!
            if (dimag(x(j)) .lt. dimag(x(j+igap))) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                goto 180
            end if
            j = j-igap
            goto 170
180         continue
        end do
        igap = igap/2
        goto 160
    end if
!
9000 continue
    call matfpe(1)
!
!     %---------------%
!     | END OF ZSORTC |
!     %---------------%
!
end subroutine
