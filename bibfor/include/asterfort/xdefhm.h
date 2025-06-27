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
    subroutine xdefhm(dimdef, dimenr, addeme, adenme, addep1,&
                      ndim, degem1, degep1, defgem, defgep, adenhy,&
                      nfh)
        integer(kind=8) :: dimenr
        integer(kind=8) :: dimdef
        integer(kind=8) :: addeme
        integer(kind=8) :: adenme
        integer(kind=8) :: addep1
        integer(kind=8) :: ndim
        real(kind=8) :: degem1(dimenr)
        real(kind=8) :: degep1(dimenr)
        real(kind=8) :: defgem(dimdef)
        real(kind=8) :: defgep(dimdef)
        integer(kind=8) :: adenhy
        integer(kind=8) :: nfh
    end subroutine xdefhm
end interface 
