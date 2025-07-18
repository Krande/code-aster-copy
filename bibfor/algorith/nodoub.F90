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
subroutine nodoub(nbl, nbb, nol, nob, typl, &
                  typb, mailla, double)
!    P. RICHARD     DATE 16/07/90
!-----------------------------------------------------------------------
!  BUT: COMPARER DEUX LISTES DE NUMEROS DE NOEUDS PAR ORDRE CROISSANT ET
    implicit none
!         DETECTER LES ELEMENTS COMMUNS
!         ARRET EN CAS D'INTERSECTION NON VIDE
!-----------------------------------------------------------------------
!
! NBL      /I/: NOMBRE DE NOEUDS INTERFACE LIBRE
! NBB      /I/: NOMBRE DE NOEUDS INTERFACE BLOQUEE
! NOL      /I/: VECTEUR DES NUMEROS DES NOEUDS LIBRES
! NOB      /I/: VECTEUR DES NUMEROS DES NOEUDS BLOQUES
! MAILLAGE /I/: NOM DU UTILISATEUR DU MAILLAGE
!
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    aster_logical :: double
    integer(kind=8) :: i, jf, lcou, lp, nbb, nbl
    character(len=24) :: valk(3)
    character(len=8) :: nomnoe, mailla, typl, typb
    integer(kind=8) :: nol(nbl), nob(nbb)
!-----------------------------------------------------------------------
!
!
    if (nbl .eq. 0 .or. nbb .eq. 0) goto 999
!
    double = .false.
    jf = 1
    do i = 1, nbl
        jf = jf-1
        lcou = nol(i)
        lp = 0
20      continue
        if (lp .lt. lcou .and. jf .lt. nbb) then
            jf = jf+1
            lp = nob(jf)
            if (lp .eq. lcou) then
                double = .true.
                nomnoe = int_to_char8(lp)
                valk(1) = nomnoe
                valk(2) = typl
                valk(3) = typb
                call utmess('E', 'ALGORITH13_69', nk=3, valk=valk)
            end if
!
            goto 20
!
        else
            goto 10
!
        end if
!
10      continue
    end do
!
    goto 999
!
999 continue
end subroutine
