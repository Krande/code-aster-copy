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
    subroutine smcarc(nb_hist, nbPhase, ftrc, trc, &
                      coef, fmod, &
                      metaSteelPara, &
                      tempCurr, tPoint, deltaTime, &
                      metaPrev, metaCurr)
        use Metallurgy_type
        integer, intent(in) :: nb_hist, nbPhase
        real(kind=8), intent(inout) :: ftrc((3*nb_hist), 3), trc((3*nb_hist), 5)
        real(kind=8), intent(in)  :: coef(*), fmod(*)
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: tempCurr, tPoint, deltaTime
        real(kind=8), intent(in) :: metaPrev(:)
        real(kind=8), intent(out) :: metaCurr(:)
    end subroutine smcarc
end interface
