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
    subroutine lc0050(BEHinteg,&
                      fami   , kpg   , ksp   , ndim  , typmod,&
                      imate  , compor, carcri, instam, instap,&
                      neps   , epsm  , deps  , nsig  , sigm  ,&
                      nvi    , vim   , option, angmas,&
                      stress , statev, dsidep, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, ndim
        character(len=8), intent(in) :: typmod(*)
        integer(kind=8), intent(in) :: imate
        character(len=16), intent(in) :: compor(*)
        real(kind=8), intent(in) :: carcri(*)
        real(kind=8), intent(in) :: instam, instap
        integer(kind=8), intent(in) :: neps, nsig, nvi
        real(kind=8), intent(in) :: epsm(6), deps(6)
        real(kind=8), intent(in) :: sigm(6)
        real(kind=8), intent(in) :: vim(*)
        character(len=16), intent(in) :: option
        real(kind=8), intent(in) :: angmas(*)
        real(kind=8), intent(out) :: stress(6)
        real(kind=8), intent(out) :: statev(nvi)
        real(kind=8), intent(out) :: dsidep(6, 6)
        integer(kind=8), intent(out) :: codret
    end subroutine lc0050
end interface
