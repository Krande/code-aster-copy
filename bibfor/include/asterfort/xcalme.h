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
    subroutine xcalme(ds_thm,&
                      option, ndim, dimenr,&
                      dimcon, addeme, adcome, congep,&
                      dsde, deps, angl_naut)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimenr
        character(len=16) :: option
        integer(kind=8) :: ndim
        integer(kind=8) :: addeme
        integer(kind=8) :: adcome
        real(kind=8) :: congep(dimcon)
        real(kind=8) :: dsde(dimcon, dimenr)
        real(kind=8) :: deps(6)
        real(kind=8), intent(in) :: angl_naut(3)
    end subroutine xcalme
end interface 
