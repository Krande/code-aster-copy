! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=C1509

subroutine lc2038(neps, nsig, option, sigp, vip, ndsde, dsidep)

    use Behaviour_type
    implicit none

#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"

    integer, intent(in) :: neps
    integer, intent(in) :: nsig
    character(len=16), intent(in) :: option
    real(kind=8), intent(out):: sigp(nsig)
    real(kind=8), intent(out):: vip(1)
    integer, intent(in) :: ndsde
    real(kind=8),      intent(out):: dsidep(merge(nsig,6,nsig*neps.eq.ndsde), merge(neps,6,nsig*neps.eq.ndsde))
! --------------------------------------------------------------------------------------------------
! Behaviour SANS
! --------------------------------------------------------------------------------------------------
    aster_logical:: lMatr, lSigm, lVari
! --------------------------------------------------------------------------------------------------

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lSigm) sigp = 0
    if (lVari) vip(1) = 0
    if (lMatr) dsidep = 0

end subroutine
