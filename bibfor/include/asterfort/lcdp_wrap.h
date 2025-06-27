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

!
!

#include "asterf_types.h"

interface
    subroutine lcdp_wrap(fami, kpg, ksp, ndim, imate, &
                        crit, instam, instap, neps, epsm,&
                        deps, vim, option, sigm, sigp, vip,& 
                        typmod, dsidep, codret)

        use lcdp_module, only: dp_material

        integer(kind=8)      :: imate, ndim, kpg, ksp, codret, neps
        real(kind=8) :: instam,instap
        real(kind=8) :: crit(*)
        real(kind=8) :: epsm(neps), deps(neps)
        real(kind=8) :: sigp(neps), sigm(neps)
        real(kind=8) :: vim(*), vip(*)
        real(kind=8) :: dsidep(neps,neps)
        character(len=16) :: option
        character(len=*) :: fami
        character(len=8) :: typmod(*)
    end subroutine 
end interface
