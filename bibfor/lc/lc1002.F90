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

subroutine lc1002(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, &
                  epsm, deps, nsig, sigm, vim, &
                  option, sigp, vip, typmod, ndsde, &
                  dsidep, codret)
!
    implicit none
!
#include "asterfort/lcpivm.h"
!
! aslint: disable=W1504,W0104
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in) :: carcri(*)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    real(kind=8), intent(in) :: vim(*)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(*)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(ndsde)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - Special SIMO_MIEHE
!
! 'VMIS_ISOT_LINE', 'VMIS_ISOT_TRAC', 'VMIS_ISOT_PUIS', 'VISC_ISOT_TRAC', 'VISC_ISOT_LINE'
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_comp
!
! --------------------------------------------------------------------------------------------------
!
    rela_comp = compor(1)
    call lcpivm(fami, kpg, ksp, imate, rela_comp, &
                carcri, instam, instap, epsm, deps, &
                vim, option, sigp, vip, dsidep, &
                codret)
!
end subroutine
