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
subroutine rrldc(a, nordre, x, nves)
!    A. COMTE                                 DATE 31/07/91
!-----------------------------------------------------------------------
!  BUT:  RESOLUTION DE L'EQUATION MATRICIELLE
    implicit none
!                 A*X=B
!
!  OU A EST UNE MATRICE COMPLEXE FACTORISEE LDLT PAR TRLDC
!
!-----------------------------------------------------------------------
!
! A        /I/: MATRICE CARRE COMPLEXE TRIANGULEE LDLT
! NORDRE   /I/: DIMENSION DE LA MATRICE A
! X        /M/: MATRICE IN:SECONDS MEMBRES   OUT:SOLUTIONS
! NVEC     /I/: NOMBRE DE VECTEURS SECOND MEMBRE
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nves, nordre
    complex(kind=8) :: a(*), x(nordre, nves), r8val
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idiag, ilign1, ilign2, in, indiag, nv
!-----------------------------------------------------------------------
    ilign1 = 1
    ilign2 = nordre
!
!     RESOLUTION DESCENDANTE
    do nv = 1, nves
        do in = ilign1, ilign2-1
            r8val = -x(in, nv)
            do i = in+1, ilign2
                idiag = i*(i-1)/2+1
                x(i, nv) = x(i, nv)+r8val*dconjg(a(idiag+i-in))
            end do
        end do
    end do
!
!     RESOLUTION DIAGONALE
    do nv = 1, nves
        do in = ilign1, ilign2
            indiag = in*(in-1)/2+1
            x(in, nv) = x(in, nv)/a(indiag)
        end do
    end do
!
!     RESOLUTION REMONTANTE
    do nv = 1, nves
        do in = ilign2, ilign1+1, -1
            indiag = in*(in-1)/2+1
            r8val = -x(in, nv)
            do i = 1, in-1
                x(i, nv) = x(i, nv)+r8val*a(indiag+in-i)
            end do
        end do
    end do
!
end subroutine
