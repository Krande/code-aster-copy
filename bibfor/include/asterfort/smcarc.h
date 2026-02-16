! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
    subroutine smcarc(nbHistTRC, ftrc, trc, &
                      coef, fmod, &
                      metaSteelPara, &
                      tempCurr, tPoint, deltaTime, &
                      nbVari, metaPrev, metaCurr)
        use Metallurgy_type
        integer(kind=8), intent(in) :: nbHistTRC
        real(kind=8), intent(inout) :: ftrc((3*nbHistTRC), 3), trc((3*nbHistTRC), 5)
        real(kind=8), intent(in)  :: coef(*), fmod(*)
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: tempCurr, tPoint, deltaTime
        integer(kind=8), intent(in) :: nbVari
        real(kind=8), intent(in) :: metaPrev(nbVari)
        real(kind=8), intent(out) :: metaCurr(nbVari)
    end subroutine smcarc
end interface
