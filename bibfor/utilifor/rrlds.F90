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
subroutine rrlds(a, nmax, nordre, x, nves)
    implicit none
!A
!A    ARGUMENTS :
!A    ---------
!A -> A      : MATRICE CARREE TRIANGULEE PAR TRLDS
!A -> NMAX   : DIMENSION REELLE DE LA MATRICE A ET DU TABLEAU X
!A -> NORDRE : DIMENSION DE LA MATRICE A
!A <> X      : SECOND MEMBRE ECRASE PAR LA SOLUTION
!A -> NVES   : NOMBRE DE COLONNES DU SECOND MEMBRE
!A
!A    BUT :
!A    ---
!A    RESOLUTION DE  A * X = B ; B ETANT STOCKE DANS X A L'APPEL
    integer(kind=8) :: nves, nmax, nordre
    real(kind=8) :: a(nmax, nordre), x(nmax, nves), r8val
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ilign1, ilign2, in, nv
!-----------------------------------------------------------------------
    ilign1 = 1
    ilign2 = nordre
!
!     RESOLUTION DESCENDANTE
    do nv = 1, nves
        do in = ilign1, ilign2-1
            r8val = -x(in, nv)
            do i = in+1, ilign2
                x(i, nv) = x(i, nv)+r8val*a(i, in)
            end do
        end do
    end do
!
!     RESOLUTION DIAGONALE
    do nv = 1, nves
        do in = ilign1, ilign2
            x(in, nv) = x(in, nv)/a(in, in)
        end do
    end do
!
!     RESOLUTION REMONTANTE
    do nv = 1, nves
        do in = ilign2, ilign1+1, -1
            r8val = -x(in, nv)
            do i = 1, in-1
                x(i, nv) = x(i, nv)+r8val*a(i, in)
            end do
        end do
    end do
!
end subroutine
