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
    subroutine vefnme_cplx(option, base, model, mate, carael, &
                           comporZ, timePrev, timeCurr, nh, ligrelz, varicomz, &
                           sigmaPrev, sigmaz, strxz, deplz, vecelz)
        character(len=16), intent(in) :: option
        character(len=1), intent(in) :: base
        character(len=8), intent(in) :: model
        real(kind=8), intent(in) :: timePrev, timeCurr
        character(len=8), intent(in) :: carael
        character(len=24), intent(in) :: mate
        character(len=*), intent(in) :: ligrelz
        integer(kind=8), intent(in) :: nh
        character(len=*), intent(in) :: comporZ
        character(len=*), intent(in) :: sigmaz, sigmaPrev
        character(len=*), intent(in) :: varicomz
        character(len=*), intent(in) :: strxz
        character(len=*), intent(in) :: deplz
        character(len=*), intent(inout) :: vecelz(*)
    end subroutine vefnme_cplx
end interface
