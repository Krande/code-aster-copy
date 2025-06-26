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

subroutine lcvisc(fami, kpg, ksp, ndim, imate, &
                  lSigm, lMatr, lVari, &
                  instam, instap, deps, vim, &
                  sigp, vip, dsidep)

    use tenseur_dime_module, only: identity, voigt

    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"

    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    integer(kind=8), intent(in)          :: ndim
    integer(kind=8), intent(in)          :: imate
    aster_logical, intent(in)   :: lSigm
    aster_logical, intent(in)   :: lMatr
    aster_logical, intent(in)   :: lVari
    real(kind=8), intent(in)     :: instam
    real(kind=8), intent(in)     :: instap
    real(kind=8), intent(in)     :: deps(:)
    real(kind=8), intent(in)     :: vim(:)
    real(kind=8), intent(inout)  :: sigp(:)
    real(kind=8), intent(out)    :: vip(:)
    real(kind=8), intent(inout)  :: dsidep(:, :)
! ----------------------------------------------------------------------
    integer(kind=8)      :: nd, iok(2)
    real(kind=8) :: k, tau, dt, a, b
    real(kind=8) :: sivm(2*ndim), siv(2*ndim), sivi(2*ndim), enerElas, incrEnerDiss
    real(kind=8) :: valev(2)
    character(len=16):: nomev(2)
! ----------------------------------------------------------------------
    data nomev/'K', 'TAU'/
! ----------------------------------------------------------------------

!   Controles
    ASSERT(size(vim) .eq. 8)
    ASSERT(size(vip) .eq. 8)
    ASSERT(size(deps) .eq. 2*ndim)
    ASSERT(size(sigp) .eq. 2*ndim)
    ASSERT(size(dsidep, 1) .eq. 2*ndim)
    ASSERT(size(dsidep, 2) .eq. 2*ndim)

!   Initialisation
    nd = 2*ndim
    dt = instap-instam
    sivm = vim(1:nd)*voigt(nd)

!   Viscous parameters
  call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'VISC_ELAS', 0, ' ', [0.d0], 2, nomev, valev, iok, 2)
    k = valev(1)
    tau = valev(2)

!   Integration scheme parameters
    a = exp(-dt/tau)
    b = k*tau/dt*(1-a)

!   Viscous stress update
    siv = a*sivm+b*deps

!   Post-treatment
    sivi = 0.5d0*(sivm+siv)
    enerElas = dot_product(siv, siv)/(2*k)
    incrEnerDiss = dot_product(sivi, sivi)/(k*tau)*dt

!   Storage
    if (lSigm) then
        sigp = sigp+siv
    end if

    if (lVari) then
        vip(1:6) = 0
        vip(1:nd) = siv/voigt(nd)
        vip(7) = enerElas
        vip(8) = vim(8)+incrEnerDiss
    end if

    if (lMatr) then
        dsidep = dsidep+b*identity(nd)
    end if

end subroutine
