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
!
interface
    subroutine lc7013(fami, kpg, ksp, ndim, imate,&
                      compor, carcri, instam, instap, epsm,&
                      deps, sigm, vim, option, angmas,&
                      sigp, vip, typmod, icomp,&
                      nvi, dsidep, codret)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ksp
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: imate
        character(len=16), intent(in) :: compor(*)
        real(kind=8), intent(in) :: carcri(*)
        real(kind=8), intent(in) :: instam
        real(kind=8), intent(in) :: instap
        real(kind=8), intent(in) :: epsm(6)
        real(kind=8), intent(in) :: deps(6)
        real(kind=8), intent(in) :: sigm(6)
        real(kind=8), intent(in) :: vim(7)
        character(len=16), intent(in) :: option
        real(kind=8), intent(in) :: angmas(*)
        real(kind=8), intent(out) :: sigp(6)
        real(kind=8), intent(out) :: vip(7)
        character(len=8), intent(in) :: typmod(*)
        integer(kind=8), intent(in) :: icomp
        integer(kind=8), intent(in) :: nvi
        real(kind=8), intent(out) :: dsidep(6, 6)
        integer(kind=8), intent(out) :: codret
    end subroutine lc7013
end interface
