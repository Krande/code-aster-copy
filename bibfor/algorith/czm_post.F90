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

subroutine czm_post(ndim, lSigm, lMatr, r, mu, su, delta, dsde, sigp, dsidep)

    implicit none
#include "asterf_types.h"

    aster_logical, intent(in) :: lMatr, lSigm
    integer(kind=8), intent(in):: ndim
    real(kind=8), intent(in):: r, mu(1:ndim), su(1:ndim), delta(1:ndim), dsde(:, :)
    real(kind=8), intent(out):: sigp(:), dsidep(:, :)
! --------------------------------------------------------------------------------------------------
!  Rangement des contraintes et des matrices tangentes pour les elements d'interface
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: i
    real(kind=8):: iden(ndim, ndim), k(ndim, ndim), m(ndim, ndim)
! --------------------------------------------------------------------------------------------------

    if (lSigm) then
        sigp(1:ndim) = mu+r*(su-delta)
        sigp(ndim+1:2*ndim) = su-delta
    end if

    if (lMatr) then
        iden = 0
        forall (i=1:ndim) iden(i, i) = 1

        k = dsde(1:ndim, 1:ndim)
        m = iden-r*k
        dsidep(1:ndim, 1:ndim) = r*m
        dsidep(ndim+1:2*ndim, 1:ndim) = m
        dsidep(1:ndim, ndim+1:2*ndim) = m
        dsidep(ndim+1:2*ndim, ndim+1:2*ndim) = -k
    end if

end subroutine
