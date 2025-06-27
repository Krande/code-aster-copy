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
!
interface
    subroutine rc36fs(nbsig1, noc1, sit1, nbsig2, noc2,&
                      sit2, saltij, ns, nscy, matse,&
                      mse, sn, nommat, c, k,&
                      cara, ug)
        integer(kind=8) :: nbsig1
        integer(kind=8) :: noc1(*)
        integer(kind=8) :: sit1(*)
        integer(kind=8) :: nbsig2
        integer(kind=8) :: noc2(*)
        integer(kind=8) :: sit2(*)
        real(kind=8) :: saltij(*)
        integer(kind=8) :: ns
        integer(kind=8) :: nscy
        real(kind=8) :: matse(*)
        real(kind=8) :: mse(*)
        real(kind=8) :: sn(*)
        character(len=8) :: nommat
        real(kind=8) :: c(*)
        real(kind=8) :: k(*)
        real(kind=8) :: cara(*)
        real(kind=8) :: ug
    end subroutine rc36fs
end interface
