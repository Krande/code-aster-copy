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
! aslint: disable=W1504,W0104,C1509

subroutine lc0019(BEHinteg, fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  ndsde, dsidep, codret)

    use Behaviour_type
    implicit none

#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/hypela.h"

    type(Behaviour_Integ)        :: BEHinteg
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
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: icomp
    integer(kind=8), intent(in) :: ndsde
    real(kind=8)                 :: dsidep(merge(nsig,6,nsig*neps.eq.ndsde), merge(neps,6,nsig*neps.eq.ndsde))
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!  RELATION DE COMPORTEMENT ELAS_HYPER
! --------------------------------------------------------------------------------------------------
    aster_logical     :: lMatr, lSigm, lVari
    integer(kind=8)           :: ndimsi
    real(kind=8)      :: eps(2*ndim), sig(6), dsde(6, 6), vi(nvi)
! --------------------------------------------------------------------------------------------------
    ASSERT(nsig .ge. 2*ndim)
    ASSERT(neps .ge. 2*ndim)

    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    ndimsi = 2*ndim
    codret = 0
    sig = 0
    vi = 0
    dsde = 0
    eps = epsm(1:ndimsi)+deps(1:ndimsi)

    call hypela(fami, kpg, ksp, ndim, typmod, imate, carcri, eps, option, sig, dsde, codret)
    if (codret .ne. 0) goto 999

    if (lSigm) sigp(1:ndimsi) = sig(1:ndimsi)
    if (lVari) vip(1:nvi) = vi(1:nvi)
    if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde(1:ndimsi, 1:ndimsi)

999 continue
end subroutine
