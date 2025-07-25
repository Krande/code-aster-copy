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

subroutine nuno1(i, iligr, nunoel, n, inum21, &
                 inuno2, nlili)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
    integer(kind=8) :: i, iligr
!
!
! ---- OBJET : FONCTION INVERSE DU CHAMP .NUNO D'UNE S.D. NUME_DDL
!      ! LE CHAMP .NUNO N'EST PLUS CONSERVE DANS LA S.D. FINALE!
! ---- DESCRIPTION DES PARAMETRES
! IN  I  I     : NUMERO DU NOEUD DANS LA NUMEROTATION .NUNO
! OUT I  ILIGR : NUMERO DANS .LILI DU LIGREL DANS LEQUEL EST DECLARE
!                LE NOEUD NUMERO I
! OUT I  NUNOEL: NUMERO DU NOEUD DANS LA NUMEROTATION LOCALE DU LIGREL
!
!-----------------------------------------------------------------------
!     FONCTIONS JEVEUX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, ili, ili1, inum21, inuno2, j
    integer(kind=8) :: jli, n, nlili, nunoel
!-----------------------------------------------------------------------
    ASSERT((i .gt. 0) .and. (i .le. n))
    j = zi(inum21+i)
    if (j .eq. 0) then
        iligr = 0
        nunoel = 0
        goto 50
    end if
!
!---- RECHERCHE DU LIGREL (CALCUL DE ILIGR)
!
    ili = 1
    i1 = zi(inuno2+ili-1)
10  continue
    do jli = ili+1, nlili+1
        i2 = zi(inuno2+jli-1)
        ili1 = jli
        if (i2 .gt. i1) exit
    end do
!
    if ((j .ge. i1) .and. (j .lt. i2)) then
        iligr = ili1-1
    else
        ili = ili1
        i1 = i2
        goto 10
    end if
!
!---- CALCUL DE NUNOEL
!
    nunoel = j-zi(inuno2+iligr-1)+1
50  continue
end subroutine
