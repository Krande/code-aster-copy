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
    subroutine epthmc(fami      , nno      , ndim  , nbsig, npg    ,&
                      shape_func, angl_naut, time , j_mater,&
                      option    , epsi_varc)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: nno
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: nbsig
        integer(kind=8), intent(in) :: npg
        real(kind=8), intent(in) :: shape_func(1)
        real(kind=8), intent(in) :: angl_naut(3)
        real(kind=8), intent(in) :: time
        integer(kind=8), intent(in) :: j_mater
        character(len=16), intent(in) :: option
        real(kind=8), intent(out) :: epsi_varc(1)
    end subroutine epthmc
end interface
