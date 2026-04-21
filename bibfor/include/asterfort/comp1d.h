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
    subroutine comp1d(BEHinteg, &
                      option, sigx, &
                      epsx, depx, vim, vip, &
                      sigxp, etan, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        character(len=16) :: option
        real(kind=8) :: sigx
        real(kind=8) :: epsx
        real(kind=8) :: depx
        real(kind=8) :: vim(*)
        real(kind=8) :: vip(*)
        real(kind=8) :: sigxp
        real(kind=8) :: etan
        integer(kind=8) :: codret
    end subroutine comp1d
end interface
