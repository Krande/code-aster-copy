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
interface
    subroutine medom2(modele, mate, mateco, cara, kcha, ncha,&
                      result, nuord, nbordr, base,&
                      npass, ligrel)
        character(len=8) :: modele
        character(len=24) :: mate, mateco
        character(len=8) :: cara
        character(len=19) :: kcha
        integer(kind=8) :: ncha
        character(len=8) :: result
        integer(kind=8) :: nuord
        integer(kind=8) :: nbordr
        character(len=1) :: base
        integer(kind=8) :: npass
        character(len=24) :: ligrel
    end subroutine medom2
end interface
