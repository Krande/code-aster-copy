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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nzgdzi(fami, kpg, ksp, ndim, imat,&
                      compor, carcri, instam, instap, fm,&
                      df, sigm, vim, option, sigp,&
                      vip, dsigdf, iret)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ksp
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: imat
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: instam
        real(kind=8), intent(in) :: instap
        real(kind=8), intent(in) :: fm(3, 3)
        real(kind=8), intent(in) :: df(3, 3)
        real(kind=8), intent(in) :: sigm(*)
        real(kind=8), intent(in) :: vim(*)
        character(len=16), intent(in) :: option
        real(kind=8), intent(out) :: sigp(*)
        real(kind=8), intent(out) :: vip(*)
        real(kind=8), intent(out) :: dsigdf(6, 3, 3)
        integer(kind=8), intent(out) :: iret
    end subroutine nzgdzi
end interface
