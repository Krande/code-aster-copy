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
subroutine ordr8(tab, nb, iord)
!    P. RICHARD     DATE //
!-----------------------------------------------------------------------
!  BUT:  TROUVER L'ORDRE CROISSANT D'UNE TABLE DE VALEUR R8
    implicit none
!     PAS DE MODIFICATION DE L'ORDRE D'ENTREE MAIS DETERMINATION DE
!     POINTEUR D'ORDRE
!
!-----------------------------------------------------------------------
!
! TAB      /I/: TABLEAU A ORDONNER
! NB       /I/: TAILLAE DU TABLEAU A ORDONNER
! IORD     /O/: TABLE DES POINTEURS D'ORDRE
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nb, iord(nb)
    real(kind=8) :: tab(nb)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iormin, itemp, j
    real(kind=8) :: vmin
!-----------------------------------------------------------------------
    do i = 1, nb
        iord(i) = i
    end do
!
    do i = 1, nb-1
        vmin = tab(iord(i))
        iormin = i
        do j = i+1, nb
            if (tab(iord(j)) .lt. vmin) then
                vmin = tab(iord(j))
                iormin = j
            end if
        end do
        itemp = iord(i)
        iord(i) = iord(iormin)
        iord(iormin) = itemp
    end do
!
end subroutine
