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
    subroutine smcarc(nb_hist, nb_phase, ftrc, trc, &
                      coef, fmod, &
                      metaSteelPara, &
                      temp_curr, temp_incr, time_incr, &
                      vari_prev, vari_curr)
        use Metallurgy_type
        integer(kind=8), intent(in) :: nb_hist, nb_phase
        real(kind=8), intent(inout) :: ftrc((3*nb_hist), 3), trc((3*nb_hist), 5)
        real(kind=8), intent(in)  :: coef(*), fmod(*)
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: temp_curr, temp_incr, time_incr
        real(kind=8), intent(in) :: vari_prev(:)
        real(kind=8), intent(out) :: vari_curr(:)
    end subroutine smcarc
end interface
