! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
    subroutine metaSteelGrainSize(metaSteelPara, &
                                  temp, time_incr1, time_incr2, &
                                  zaustenite, coef_phase, &
                                  d_prev, d_curr)
        use Metallurgy_type
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: d_prev, temp, time_incr1, time_incr2
        real(kind=8), intent(in) :: zaustenite, coef_phase
        real(kind=8), intent(out) :: d_curr
    end subroutine metaSteelGrainSize
end interface
