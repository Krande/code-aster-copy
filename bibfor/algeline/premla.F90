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
subroutine premla(neq, diag, col, lt, nrl, &
                  rl, deb, vois, suit, ier)
! person_in_charge: olivier.boiteau at edf.fr
!
    implicit none
#include "asterfort/calajt.h"
#include "asterfort/infniv.h"
    integer(kind=8) :: neq, diag(0:neq), col(*), deb(neq)
    integer(kind=8) :: vois(*), suit(*), lt, nrl, rl(4, nrl), ier
!     VARIABLES LOCALES
    integer(kind=8) :: i, j, j1, j2, k, illist, ifm, niv, lbd2
!****************************************************************
!-----RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, niv)
    ier = 0
    if (nrl .ne. 0) then
!
!     AVEC RELATION LINEAIRE
!     CALCUL DES LISTES DE NOEUDS A AJOUTER
        do i = 1, neq
            deb(i) = 0
        end do
        illist = 0
        do k = 1, nrl
            lbd2 = rl(2, k)
            j1 = diag(lbd2-1)+2
            j2 = diag(lbd2)-1
            do j = j2, j1, -1
!     ON AJOUTE COL(J1),..., COL(J-1) AUX VOISINS DE COL(J)
                call calajt(j1, j, diag, col, neq, &
                            illist, deb, vois, suit, lt, &
                            ier)
                if (ier .gt. 0) goto 999
            end do
!
        end do
    end if
999 continue
!
!
end subroutine
