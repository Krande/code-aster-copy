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

subroutine lc3053(BEHinteg, fami, kpg, ksp, ndim, imate, &
                  compor, crit, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  ndsde, dsidep, codret)
! aslint: disable=W1504,W0104

    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lckimp.h"

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
    real(kind=8)                 :: dsidep(nsig, neps)
    integer(kind=8), intent(out):: codret
! --------------------------------------------------------------------------------------------------
!  RELATION DE COMPORTEMENT ENDO_CARRE pour modelisation GVNO
! --------------------------------------------------------------------------------------------------
    aster_logical     :: lMatr, lSigm, lVari
    integer(kind=8)           :: ndimsi
    real(kind=8)      :: eps(2*ndim), phi, sig(2*ndim), forc_endo, vi(nvi)
    real(kind=8)      :: dsde_1(2*ndim, 2*ndim), dsde_2(2*ndim), dsde_3
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nsig .ge. 2*ndim+1)
    ASSERT(neps .ge. 2*ndim+1)
    ASSERT(ndsde .ge. neps*nsig)

    ndimsi = 2*ndim
    codret = 0
    sig = 0
    forc_endo = 0
    vi = 0
    dsde_1 = 0
    dsde_2 = 0
    dsde_3 = 0

    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

    eps = epsm(1:ndimsi)+deps(1:ndimsi)
    phi = epsm(ndimsi+1)+deps(ndimsi+1)

    call lckimp(ndim, typmod(1), option, imate, eps, &
                phi, vim, sig, forc_endo, vi, dsde_1, dsde_2, dsde_3)

    if (lSigm) then
        sigp(1:ndimsi) = sig
        sigp(ndimsi+1) = forc_endo
    end if

    if (lVari) vip(1:nvi) = vi

    if (lMatr) then
        dsidep(1:ndimsi, 1:ndimsi) = dsde_1
        dsidep(ndimsi+1, 1:ndimsi) = dsde_2
        dsidep(1:ndimsi, ndimsi+1) = dsde_2
        dsidep(ndimsi+1, ndimsi+1) = dsde_3
    end if

end subroutine
