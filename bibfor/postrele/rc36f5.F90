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
subroutine rc36f5(nbp12, nbp23, nbp13, nbsigr, nbsg1, &
                  nbsg2, nbsg3, saltij)
    implicit none
    integer(kind=8) :: nbp12, nbp23, nbp13, nbsigr, nbsg1, nbsg2, nbsg3
    real(kind=8) :: saltij(*)
!
!     SI IL N'EXISTE PAS DE SITUATION DE PASSAGE ENTRE 2 GROUPES,
!     ON MET LES TERMES CROISES DE SALT A ZERO
!
!     ------------------------------------------------------------------
    integer(kind=8) :: i1, i2, isl
!     ------------------------------------------------------------------
!
    if (nbp12 .eq. 0) then
!           BLOC 1_2
        do i1 = 1, nbsg1
            isl = 4*(i1-1)*nbsigr+4*nbsg1
            do i2 = 1, nbsg2
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
!           BLOC 2_1
        do i1 = 1, nbsg2
            isl = 4*nbsigr*nbsg1+4*(i1-1)*nbsigr
            do i2 = 1, nbsg1
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
    end if
!
    if (nbp23 .eq. 0) then
!           BLOC 2_3
        do i1 = 1, nbsg2
            isl = 4*nbsigr*nbsg1+4*(i1-1)*nbsigr+4*(nbsg1+nbsg2)
            do i2 = 1, nbsg3
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
!           BLOC 3_2
        do i1 = 1, nbsg3
            isl = 4*nbsigr*(nbsg1+nbsg2)+4*(i1-1)*nbsigr+4*nbsg1
            do i2 = 1, nbsg2
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
    end if
!
    if (nbp13 .eq. 0) then
!           BLOC 1_3
        do i1 = 1, nbsg1
            isl = 4*(i1-1)*nbsigr+4*(nbsg1+nbsg2)
            do i2 = 1, nbsg3
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
!           BLOC 3_1
        do i1 = 1, nbsg3
            isl = 4*nbsigr*(nbsg1+nbsg2)+4*nbsigr*(i1-1)
            do i2 = 1, nbsg1
                saltij(isl+4*(i2-1)+1) = 0.d0
                saltij(isl+4*(i2-1)+2) = 0.d0
                saltij(isl+4*(i2-1)+3) = 0.d0
                saltij(isl+4*(i2-1)+4) = 0.d0
            end do
        end do
    end if
!
end subroutine
