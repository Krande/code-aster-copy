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
    subroutine lcejmr(BEHinteg, &
                      fami, kpg, ksp, ndim, mate, &
                      option, epsm, deps, sigma, &
                      dsidep, vim, vip, typmod, &
                      instam, instap)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer(kind=8) :: ndim
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: mate
        character(len=16) :: option
        real(kind=8) :: epsm(6)
        real(kind=8) :: deps(6)
        real(kind=8) :: sigma(6)
        real(kind=8) :: dsidep(6, 6)
        real(kind=8) :: vim(*)
        real(kind=8) :: vip(*)
        character(len=8) :: typmod(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
    end subroutine lcejmr
end interface
