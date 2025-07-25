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
subroutine trlds(a, nmax, nordre, ierr)
    implicit none
!A
!A    ARGUMENTS :
!A    ---------
!A <> A      : MATRICE CARREE
!A -> NMAX   : DIMENSION REELLE DE LA MATRICE A
!A -> NORDRE : DIMENSION DE LA MATRICE A
!A <- IERR   : NUMERO DE LA LIGNE POUR LAQUELLE UN PIVOT NUL EST OBTENU
!A             0 SI PAS DE PIVOT NUL
!A
!A    BUT :
!A    ---
!A    TRIANGULATION EN PLACE DE LA MATRICE CARREE A
#include "asterfort/iunifi.h"
    integer(kind=8) :: nmax, nordre
    real(kind=8) :: a(nmax, nordre), r8val
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibm, ierr, ifm, in, jn
!-----------------------------------------------------------------------
    ierr = 0
    do in = 1, nordre
        if (in .eq. 1) goto 50
!
!        UTILISATION  DES  LIGNES  (1) A (IN-1)
        do jn = 1, in-1
!
            if (jn .eq. 1) goto 36
            ibm = jn-1
!
            r8val = a(jn, in)
            do i = 1, ibm
                r8val = r8val-a(jn, i)*a(i, in)*a(i, i)
            end do
            a(jn, in) = r8val
!
            r8val = a(in, jn)
            do i = 1, ibm
                r8val = r8val-a(in, i)*a(i, jn)*a(i, i)
            end do
            a(in, jn) = r8val
!
36          continue
            a(jn, in) = a(jn, in)/a(jn, jn)
            a(in, jn) = a(in, jn)/a(jn, jn)
        end do
!
50      continue
!
!        UTILISATION  DE LA LIGNE IN ( CALCUL DU TERME PIVOT)
        ibm = in-1
!
        r8val = a(in, in)
        do i = 1, ibm
            r8val = r8val-a(in, i)*a(i, in)*a(i, i)
        end do
        a(in, in) = r8val
        if (r8val .eq. 0.d00) then
            ierr = in
            ifm = iunifi('MESSAGE')
            write (ifm, *) ' TRLDS : PIVOT NUL A LA LIGNE ', in
            goto 999
        end if
!
    end do
!
999 continue
end subroutine
