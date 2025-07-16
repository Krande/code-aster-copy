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
subroutine trnuli(itab, nblig, nbcol, icol, nures)
    implicit none
#include "asterf_types.h"
!
!
!***********************************************************************
!    P. RICHARD     DATE 20/01/92
!-----------------------------------------------------------------------
!  BUT:  < TROUVER NUMERO DE  LIGNE DANS UN TABLEAU >
!
!   A PARTIR DES VALEURS DES COLONNES
!
!-----------------------------------------------------------------------
!
! NOM----- / /:
!
! ITAB     /I/: TABLEAU D'ENTIER
! NBLIG    /I/: NOMBRE DE LIGNES
! NBCOL    /I/: NOMBRE DE COLONNES
! ICOL     /I/: VALEURS DES COLONNES A TROUVER
! NURES    /I/: NUMERO DE LA LIGNE CHERCHEE
!
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, nbcol, nblig, nures
    integer(kind=8) :: itab(nblig, nbcol), icol(nbcol)
    aster_logical :: ok
!-----------------------------------------------------------------------
    i = 0
    nures = 0
!
10  continue
    i = i+1
!
    ok = .true.
!
    do j = 1, nbcol
        if (itab(i, j) .ne. icol(j)) ok = .false.
    end do
!
    if (ok) then
        nures = i
        goto 999
    else
        if (i .lt. nblig) then
            goto 10
        else
            goto 999
        end if
    end if
!
!
999 continue
end subroutine
