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
! aslint: disable=W1504
!
subroutine lc0137(BEHinteg, &
                  fami, kpg, ksp, ndim, imate, &
                  compor, mult_comp, carcri, instam, instap, &
                  neps, epsm, deps, sigm, vim, option, &
                  angmas, sigp, vip, &
                  typmod, icomp, &
                  nvi, dsidep, codret)
!
    use Behaviour_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/nmvprk.h"
#include "asterfort/plasti.h"
#include "asterfort/utlcal.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(*)
    character(len=16), intent(in) :: mult_comp
    real(kind=8), intent(in) :: carcri(*)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    real(kind=8), intent(in) :: sigm(6)
    real(kind=8), intent(in) :: vim(*)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: angmas(3)
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8), intent(out) :: vip(*)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: icomp
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! polycristal, monocristal
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHinteg       : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algo_inte
    character(len=11) :: meting
    common/meti/meting
!
! --------------------------------------------------------------------------------------------------
!
    if (compor(1) .eq. 'POLYCRISTAL') then
        call nmvprk(fami, kpg, ksp, ndim, typmod, &
                    imate, compor, carcri, instam, instap, &
                    neps, epsm, deps, sigm, nvi, vim, &
                    option, angmas, sigp, vip, dsidep, &
                    codret, mult_comp)
    elseif (compor(1) .eq. 'MONOCRISTAL') then
        call utlcal('VALE_NOM', algo_inte, carcri(6))
        if (algo_inte(1:6) .eq. 'NEWTON') then
            meting = algo_inte(1:11)
            call plasti(BEHinteg, &
                        fami, kpg, ksp, typmod, imate, &
                        compor, carcri, instam, instap, &
                        epsm, deps, sigm, &
                        vim, option, angmas, sigp, vip, &
                        dsidep, icomp, nvi, codret, mult_comp)
        else if (algo_inte .eq. 'RUNGE_KUTTA') then
            meting = 'RUNGE_KUTTA'
            call nmvprk(fami, kpg, ksp, ndim, typmod, &
                        imate, compor, carcri, instam, instap, &
                        neps, epsm, deps, sigm, nvi, vim, &
                        option, angmas, sigp, vip, dsidep, &
                        codret, mult_comp)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
