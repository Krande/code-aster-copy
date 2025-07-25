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
! aslint: disable=W0104
!
subroutine lc2001(BEHinteg, &
                  fami, kpg, ksp, ndim, imate, &
                  neps, deps, nsig, sigm, option, &
                  angmas, sigp, vip, typmod, ndsde, &
                  dsidep, codret)
!
    use Behaviour_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/nmelas_incr.h"
#include "asterfort/rccoma.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: angmas(3)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(1)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(nsig, neps)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - Special IMPLEX
!
! 'ELAS'
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: mcmate
    integer(kind=8) :: icodre
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(neps*nsig .eq. ndsde)
    ASSERT(neps .eq. nsig)
    ASSERT(neps .ge. 2*ndim)

    codret = 0
    call rccoma(imate, 'ELAS', 1, mcmate, icodre)
    ASSERT(icodre .eq. 0)
    if (mcmate .eq. 'ELAS') then
        call nmelas_incr(BEHinteg, &
                         fami, kpg, ksp, typmod, &
                         imate, deps(1:2*ndim), sigm(1:2*ndim), option, &
                         sigp(1:2*ndim), vip, dsidep(1:2*ndim, 1:2*ndim))
    else
        ASSERT(.false.)
    end if
!
end subroutine
