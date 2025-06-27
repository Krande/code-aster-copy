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
#include "asterf_types.h"
!
interface
    subroutine nmcomp(BEHinteg, &
                      fami,  kpg,    ksp,    ndim,       typmod,       &
                      imate, compor, carcri, instam,     instap,       &
                      neps,  epsm_inp,   deps_inp,   nsig,       sigm,         &
                      vim,   option, angmas, sigp,       vip,          &
                      ndsde, dsidep, codret, mult_comp_, l_epsi_varc_, &
                      materi_ )
!
        use Behaviour_type
!
        type(Behaviour_Integ) :: BEHinteg
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: ndim
        character(len=8) :: typmod(*)
        integer(kind=8) :: imate
        character(len=16) :: compor(*)
        real(kind=8) :: carcri(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
        integer(kind=8) :: neps
        real(kind=8) :: epsm_inp(neps)
        real(kind=8) :: deps_inp(neps)
        integer(kind=8) :: nsig
        real(kind=8) :: sigm(nsig)
        real(kind=8) :: vim(*)
        character(len=16) :: option
        real(kind=8) :: angmas(*)
        real(kind=8) :: sigp(nsig)
        real(kind=8) :: vip(*)
        integer(kind=8) :: ndsde
        real(kind=8) :: dsidep(merge(nsig,6,nsig*neps.eq.ndsde),merge(neps,6,nsig*neps.eq.ndsde))
        integer(kind=8) :: codret
        character(len=16), optional, intent(in) :: mult_comp_
        aster_logical, optional, intent(in)     :: l_epsi_varc_
        character(len=8), optional, intent(in)  :: materi_
    end subroutine nmcomp
end interface
