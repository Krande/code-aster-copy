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
    subroutine lcmmjv(mult_comp, nmat, cpmono, nbfsys, irota,&
                      itbint, nsg, hsr)
        integer(kind=8) :: nsg
        integer(kind=8) :: nmat
        character(len=16) :: mult_comp
        character(len=24) :: cpmono(5*nmat+1)
        integer(kind=8) :: nbfsys
        integer(kind=8) :: irota
        integer(kind=8) :: itbint
        real(kind=8) :: hsr(nsg, nsg)
    end subroutine lcmmjv
end interface
