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

subroutine lc0152(BEHinteg, fami, kpg, ksp, ndim, imate, &
                  compor, crit, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  ndsde, dsidep, codret)
! aslint: disable=W1504,W0104,C1509

    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lccgad.h"

    type(Behaviour_Integ)        :: BEHinteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in) :: crit(*)
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
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: icomp
    integer(kind=8), intent(in) :: ndsde
    real(kind=8)                 :: dsidep(merge(nsig,6,nsig*neps.eq.ndsde), merge(neps,6,nsig*neps.eq.ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!  RELATION DE COMPORTEMENT CABLE_GAINE
! --------------------------------------------------------------------------------------------------
    aster_logical :: lMatr, lSigm, lVari
    real(kind=8)  :: mu, su, eps(2), de, dsde(2), vi(nvi)
! --------------------------------------------------------------------------------------------------
    ASSERT(neps .ge. 2)

    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lVari = L_VARI(option)

    codret = 0
    de = 0
    dsde = 0
    vi = 0

    eps = epsm(1:2)+deps(1:2)
    mu = eps(1)
    su = eps(2)

    call lccgad(BEHinteg, &
                fami, kpg, ksp, imate, option, &
                mu, su, de, dsde, vim, vi)

    if (lSigm) sigp(1) = de
    if (lMatr) dsidep(1, 1:2) = dsde
    if (lVari) vip(1:nvi) = vi

end subroutine
