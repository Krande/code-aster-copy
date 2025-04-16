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
    subroutine dilcar(option, icontm, ivarim, ideplm, ideplp,&
                      igeom, imate, imatuu, ivectu, icontp,&
                      ivarip, ichg, ichn, jcret, icarcr, iinstm, iinstp)
        character(len=16) :: option
        integer :: icontm
        integer :: ivarim
        integer :: ideplm
        integer :: ideplp
        integer :: igeom
        integer :: imate
        integer :: imatuu
        integer :: ivectu
        integer :: icontp
        integer :: ivarip
        integer :: ichg
        integer :: ichn
        integer :: jcret
        integer :: icarcr
        integer :: iinstm
        integer :: iinstp
    end subroutine dilcar
end interface
