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
    subroutine cgDiscrField(cgField, cgTheta, cgStudy, cgStat, chsdeg, chslag, v_absc, v_basf, &
     v_cesv, jcesd, jcesl, i_theta, lpain, lchin, nchin)
use calcG_type
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_Study), intent(in) :: cgStudy
        type(CalcG_stat), intent(inout) :: cgStat
        character(len=19), intent(in) :: chsdeg, chslag
        integer(kind=8), intent(in) :: jcesd, jcesl, i_theta
        real(kind=8), pointer :: v_basf(:)
        real(kind=8), pointer :: v_absc(:)
        integer(kind=8), pointer  :: v_cesv(:)
        character(len=24), intent(inout) :: lchin(*)
        character(len=8), intent(inout) :: lpain(*)
        integer(kind=8), intent(inout) :: nchin
    end subroutine cgDiscrField
end interface
