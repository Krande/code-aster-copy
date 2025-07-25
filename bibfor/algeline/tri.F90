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

subroutine tri(clef, tab, ntab, n)
!
! A_UTIL
! ----------------------------------------------------------------------
!                     TRI RAPIDE (HOARE / SEDGEWICK)
! ----------------------------------------------------------------------
! VARIABLES D'ENTREE / SORTIE
! INTEGER CLEF(N)         : VECTEUR CLEF
! INTEGER TAB(N,NTAB)     : TABLEAU A TRIER EN MEME TEMPS QUE CLEF
!                           (SI NTAB = 0, PAS PRIS EN COMPTE)
!
! VARIABLES D'ENTREE
! INTEGER NTAB            : NOMBRE DE COLONNES DE TAB
! INTEGER N               : NOMBRE DE LIGNES A TRIER
! ----------------------------------------------------------------------
!
    implicit none
!
! --- PARAMETRES
#include "asterfort/assert.h"
#include "asterfort/triins.h"
#include "asterfort/trirap.h"
    integer(kind=8) :: blocmx, npile
    parameter(blocmx=14)
    parameter(npile=59)
!
! --- VARIABLES
    integer(kind=8) :: n, ntab, clef(*), tab(n, *)
    integer(kind=8) :: pile(npile+1), g, d, gs, ds, m, ipile
!
! --- INITIALISATION
!
    if (n .le. blocmx) goto 20
!
    g = 1
    d = n
    ipile = 1
!
10  continue
!
! --- DECOUPAGE
!
    call trirap(clef, tab, ntab, n, g, &
                d, m)
!
    if ((m-g) .gt. (d-m)) then
        gs = g
        ds = m-1
        g = m+1
    else
        gs = m+1
        ds = d
        d = m-1
    end if
!
    if ((d-g) .ge. blocmx) then
!
! ----- PUSH
!
        if ((ds-gs) .ge. blocmx) then
!
            if (ipile .le. npile) then
!
                pile(ipile) = gs
                ipile = ipile+1
                pile(ipile) = ds
                ipile = ipile+1
!
            else
!
                ASSERT(.false.)
!
            end if
!
        end if
!
        goto 10
!
    else
!
! ----- POP
!
        if (ipile .gt. 2) then
!
            ipile = ipile-1
            d = pile(ipile)
            ipile = ipile-1
            g = pile(ipile)
            goto 10
!
        end if
!
    end if
!
! --- TRI PAR INSERTION
!
20  continue
!
    call triins(clef, tab, ntab, n)
!
end subroutine
