! --------------------------------------------------------------------
! Copyright (C) LINPACK
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
function ldasum(n, dx, incx)
!     SUBROUTINE LINPACK CALCULANT UNE SOMME DE VALEUR ABSOLUE.
!-----------------------------------------------------------------------
!     TAKES THE SUM OF THE ABSOLUTE VALUES.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!     MODIFIED TO CORRECT PROBLEM WITH NEGATIVE INCREMENT, 8/21/90.
!
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 2 RETURN PAR GOTO 1000,
!            REMPLACEMENT DE DABS PAR ABS,
!            IMPLICIT NONE.
! INTRINSIC FUNCTIONS
!            ABS.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
    real(kind=8) :: ldasum
!
    real(kind=8) :: dx(1), dtemp
    integer(kind=8) :: i, incx, ix, m, mp1, n
!
    ldasum = 0.0d0
    dtemp = 0.0d0
    if (n .le. 0) goto 1000
    if (incx .eq. 1) goto 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
    ix = 1
    if (incx .lt. 0) ix = (-n+1)*incx+1
    do i = 1, n
        dtemp = dtemp+abs(dx(ix))
        ix = ix+incx
    end do
    ldasum = dtemp
    goto 1000
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
! DUE TO CRP_11
20  continue
    m = mod(n, 6)
    if (m .eq. 0) goto 40
    do i = 1, m
        dtemp = dtemp+abs(dx(i))
    end do
    if (n .lt. 6) goto 60
! DUE TO CRP_11
40  continue
    mp1 = m+1
    do i = mp1, n, 6
        dtemp = dtemp+abs( &
                dx(i))+abs(dx(i+1))+abs(dx(i+2))+abs(dx(i+3))+abs(dx(i+4))+abs(&
                &dx(i+5) &
                )
    end do
! DUE TO CRP_11
60  continue
    ldasum = dtemp
!
1000 continue
end function
