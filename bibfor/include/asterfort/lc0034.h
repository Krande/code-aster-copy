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
    subroutine lc0034(fami, kpg, ksp, imate,&
                      carcri, epsm,&
                      deps, sigm, vim, option, angmas,&
                      sigp, vip, typmod, icomp,&
                      dsidep, codret)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ksp
        integer(kind=8), intent(in) :: imate
        real(kind=8) :: carcri(*)
        real(kind=8), intent(in) :: epsm(*)
        real(kind=8), intent(in) :: deps(*)
        real(kind=8), intent(in) :: sigm(6)
        real(kind=8) :: vim(50)
        character(len=16), intent(in) :: option
        real(kind=8), intent(in) :: angmas(3)
        real(kind=8), intent(out) :: sigp(6)
        real(kind=8) :: vip(50)
        character(len=8), intent(in) :: typmod(*)
        integer(kind=8), intent(in) :: icomp
        real(kind=8), intent(out) :: dsidep(6, 6)
        integer(kind=8), intent(out) :: codret
    end subroutine lc0034
end interface
