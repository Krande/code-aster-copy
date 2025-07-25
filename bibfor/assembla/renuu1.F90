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
subroutine renuu1(coin, longi, ordo, longo, nbco, &
                  newn)
    implicit none
!
!     ARGUMENTS:
!     ----------
    integer(kind=8) :: longi, longo, coin(*), ordo(*), nbco(*), newn(*)
! ----------------------------------------------------------------------
!     BUT:     ORDONNER LA LISTE DE NOEUDS COIN(*) DE LONGUEUR LONGI
!           DANS LA LISTE ORDO(*) :
!              UN NOEUD I EST PLACE AVANT UN NOEUD J SI SON NOMBRE DE
!           CONNECTIVITE EST INFERIEUR A CELUI DE J.
!
!              LES NOEUDS DEJA RENUMEROTES NE SONT PAS TRAITES
!           (CE SONT LES NOEUDS DONT LE .NEWN EST DEJA >0)
!           LA LISTE DES NOEUDS ORDONNES PEUT DONC ETRE DE LONGUEUR
!           MOINDRE QUE COIN : LONGO <= LONGI)
!
!     IN : LONGI : LONGUEUR UTILE DE COIN
!          COIN(*) : LISTE DES NOEUDS A ORDONNER.
!          NBCO(I) : NOMBRE DE NOEUDS CONNECTES AU NOEUD I.
!          NEWN(I) : NOUVEAU NUMERO DU NOEUD I (0 SI PAS RENUMEROTE).
!
!     OUT: LONGO : LONGUEUR UTILE DE ORDO
!          ORDO(*) : LISTE DES NOEUDS REORDONNES.
!
! ----------------------------------------------------------------------
! DEB-------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
!-----------------------------------------------------------------------
    longo = 0
    do i = 1, longi
        if (newn(coin(i)) .gt. 0) goto 1
        longo = longo+1
!
!        -- ON PLACE LE 1ER EN 1ER :
        if (longo .eq. 1) then
            ordo(1) = coin(i)
            goto 1
        end if
!
        do j = 1, longo-1
!           -- ON INSERE COIN(I) :
            if (nbco(coin(i)) .le. nbco(ordo(j))) then
                do k = longo-1, j, -1
                    ordo(k+1) = ordo(k)
                end do
                ordo(j) = coin(i)
                goto 1
            end if
        end do
        ordo(longo) = coin(i)
1       continue
    end do
!
end subroutine
