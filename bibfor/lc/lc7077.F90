! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine lc7077(BEHinteg,&
                  fami, kpg, ksp, ndim, imate,&
                  compor, carcri, instam, instap, neps, epsm,&
                  deps, nsig, sigm, nvi, vim, option, angmas,&
                  sigp, vip, typmod, icomp,&
                  ndsde, dsidep, codret)
!
use Behaviour_type
use Behaviour_module, only: behaviourOption
use czm_elas_module,  only: CONSTITUTIVE_LAW, Init, Integrate 
!
implicit none
#include "asterfort/assert.h"
!
!
! aslint: disable=W1504,W0104
!
    type(Behaviour_Integ), intent(inout) :: BEHinteg
    character(len=*), intent(in)  :: fami
    integer, intent(in)           :: kpg
    integer, intent(in)           :: ksp
    integer, intent(in)           :: ndim
    integer, intent(in)           :: imate
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in)      :: carcri(*)
    real(kind=8), intent(in)      :: instam
    real(kind=8), intent(in)      :: instap
    integer, intent(in)           :: neps
    real(kind=8), intent(in)      :: epsm(neps)
    real(kind=8), intent(in)      :: deps(neps)
    integer, intent(in)           :: nsig
    real(kind=8), intent(in)      :: sigm(nsig)
    integer, intent(in)           :: nvi
    real(kind=8), intent(in)      :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in)      :: angmas(3)
    real(kind=8), intent(out)     :: sigp(nsig)
    real(kind=8), intent(out)     :: vip(nvi)
    character(len=8), intent(in)  :: typmod(*)
    integer, intent(in)           :: icomp
    integer, intent(in)           :: ndsde
    real(kind=8), intent(out)     :: dsidep(nint(sqrt(ndsde*1.d0)),nint(sqrt(ndsde*1.d0)))
    integer, intent(out)          :: codret
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! CZM_ELAS
!
! --------------------------------------------------------------------------------------------------
    aster_logical :: lMatr, lVect, lSigm, lVari
    real(kind=8)  :: t(ndim), su(ndim), delta(ndim), dphi_delta(ndim,ndim), vi(nvi)
    type(CONSTITUTIVE_LAW):: cl
! --------------------------------------------------------------------------------------------------
    ASSERT (neps  .ge. ndim)
    ASSERT (nsig  .ge. ndim)
    ASSERT (ndsde .ge. ndim*ndim)
    ASSERT (nint(sqrt(ndsde*1.d0))*nint(sqrt(ndsde*1.d0)) .eq. ndsde)
! --------------------------------------------------------------------------------------------------
    
    t  = epsm(1:ndim)
    su = deps(1:ndim)
    
    call behaviourOption(option, compor,lMatr , lVect ,lVari , lSigm)

    cl = Init(ndim, fami, kpg, ksp, imate, t, su)
    
    call Integrate(cl, delta, dphi_delta, vi)

    codret = cl%exception
    if (codret.ne.0) goto 999

    if (lSigm) then 
        sigp = 0
        sigp(1:ndim) = delta
    end if
    
    if (lVari) vip = vi
    
    if (lMatr) then
        dsidep = 0
        dsidep(1:ndim,1:ndim) = dphi_delta
    end if

    BEHinteg%elga%r = cl%r

999 continue                      
end subroutine
 
