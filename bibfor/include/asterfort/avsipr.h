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
    subroutine avsipr(nbordr, vwork, tdisp, kwork, sommw,&
                      tspaq, i, jvsipr, jvepsn)
        integer(kind=8) :: tdisp
        integer(kind=8) :: nbordr
        real(kind=8) :: vwork(tdisp)
        integer(kind=8) :: kwork
        integer(kind=8) :: sommw
        integer(kind=8) :: tspaq
        integer(kind=8) :: i
        integer(kind=8) :: jvsipr
        integer(kind=8) :: jvepsn
    end subroutine avsipr
end interface
