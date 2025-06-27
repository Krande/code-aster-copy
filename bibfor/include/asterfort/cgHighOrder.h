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
    subroutine cgHighOrder(cgField, cgTheta, cgStudy, cgStat, chscer, chseli, v_cer, v_eli, &
     jcesd2, jcesl2, chcer, cheli)
use calcG_type
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_Study), intent(in) :: cgStudy
        type(CalcG_stat), intent(inout) :: cgStat
        character(len=19), intent(in) :: chscer, chseli
        integer(kind=8), intent(in) :: jcesd2, jcesl2
        real(kind=8), pointer :: v_cer(:)
        real(kind=8), pointer :: v_eli(:)
        character(len=19), intent(inout) :: chcer
        character(len=19), intent(inout) :: cheli
    end subroutine cgHighOrder
end interface
