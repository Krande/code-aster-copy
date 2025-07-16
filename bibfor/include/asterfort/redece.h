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
! aslint: disable=C1505

#include "asterf_types.h"
!
interface
    subroutine redece(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  l_epsi_varc, imate, materi, compor, mult_comp, &
                  carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  angmas, numlc, sigp, vip, &
                  ndsde, dsidep, codret)
        use Behaviour_type
    type(Behaviour_Integ) :: BEHinteg
    character(len=*) :: fami
    integer(kind=8) :: imate, ndim, kpg, ksp, numlc
    integer(kind=8) :: neps, nsig, ndsde
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8), intent(in):: nvi
    real(kind=8) :: carcri(*), angmas(*)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(neps), deps(neps)
    real(kind=8) :: sigm(nsig), sigp(nsig)
    real(kind=8) :: vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                              merge(neps, 6, nsig*neps .eq. ndsde))
    character(len=8)  :: typmod(*)
    character(len=16) :: compor(*), option
    character(len=8), intent(in) :: materi
    character(len=16), intent(in) :: mult_comp
    integer(kind=8),intent(out):: codret
    end subroutine redece
end interface
