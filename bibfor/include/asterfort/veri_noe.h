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
    subroutine veri_noe(mailla, dmax_cable, lnuma, liproj,&
                        nbmaok, x3dca, iproj, noe, numail )
        character(len=8) :: mailla
        real(kind=8) :: dmax_cable
        integer(kind=8) :: lnuma(*)
        integer(kind=8) :: liproj(*)
        integer(kind=8) :: nbmaok
        real(kind=8) :: x3dca(3)
        integer(kind=8) :: iproj
        integer(kind=8) :: noe
        integer(kind=8) :: numail
    end subroutine veri_noe
end interface
