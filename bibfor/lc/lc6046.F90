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

subroutine lc6046(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, &
                  epsm, deps, nsig, sigm, nvi, vim, &
                  option, angmas, sigp, vip, &
                  typmod, icomp, ndsde, &
                  dsidep, codret)
!
    implicit none
!
#include "asterfort/lcesgv.h"
#include "asterfort/lcquma.h"
#include "asterfort/lcquga.h"
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
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: angmas(*)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: icomp
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(ndsde)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - Special GRADVARI
!
! ENDO_SCALAIRE
!
! --------------------------------------------------------------------------------------------------
!
    call lcesgv(fami, kpg, ksp, ndim, neps, typmod, option, imate, lcquma, lcquga, &
                epsm, deps, vim, nint(carcri(1)), carcri(3), sigp, &
                vip, dsidep, codret)

!
end subroutine
