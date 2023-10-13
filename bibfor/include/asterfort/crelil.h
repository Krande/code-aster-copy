! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
    subroutine crelil(kstop, nbmat, vlimat, lili, base,&
                      nomma, pref, gd, mailla, nec,&
                      ncmp, ilimo, nlili, nbelm, nume_)
        character(len=1) :: kstop
        integer :: nbmat
        character(len=*) :: vlimat(nbmat)
        character(len=*) :: lili
        character(len=1) :: base
        character(len=*) :: nomma
        character(len=*) :: pref
        integer :: gd
        character(len=*) :: mailla
        integer :: nec
        integer :: ncmp
        integer :: ilimo
        integer :: nlili
        integer :: nbelm
        character(len=14), optional :: nume_
    end subroutine crelil
end interface
