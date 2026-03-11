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
! aslint: disable=W1504
!
subroutine lc0031(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, &
                  epsm, deps, sigm, nvi, vim, option, &
                  angmas, sigp, vip, typmod, &
                  dsidep, codret)
!
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmveei.h"
#include "asterfort/nmvprk.h"
#include "asterfort/utlcal.h"
!
    integer(kind=8) :: imate, ndim, kpg, ksp, codret, nvi, neps
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: angmas(*), instam, instap
    real(kind=8) :: epsm(6), deps(6), sigm(6), sigp(6), vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(6, 6)
    character(len=16) :: option
    character(len=8) :: typmod(*)
    character(len=*) :: fami
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! VENDOCHAB / VISC_ENDO_LEMA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algo_inte
!
! --------------------------------------------------------------------------------------------------
!
    call utlcal('VALE_NOM', algo_inte, carcri(6))
    if (algo_inte .eq. 'RUNGE_KUTTA') then
        call nmvprk(fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    neps, epsm, deps, sigm, nvi, vim, &
                    option, angmas, sigp, vip, dsidep, &
                    codret)
    else
        call nmveei(fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, dsidep, codret)
    end if
!
end subroutine
