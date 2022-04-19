! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
interface
    subroutine compor_3d_fibre(fami,    kpg,    ksp,  option, sigx,   &
                               instam,  instap, crit, icdmat, materi, &
                               pcompor, epsx,   depx, angmas, vim,    &
                               vip,     sigxp,  etan, codret)
        character(len=*)    :: fami
        integer             :: kpg, ksp
        character(len=16)   :: option
        real(kind=8)        :: sigx, instam, instap, crit(*)
        integer             :: icdmat
        character(len=8)    :: materi
        character(len=16)   :: pcompor(*)
        real(kind=8)        :: epsx, depx
        real(kind=8)        :: angmas(3)
        real(kind=8)        :: vim(*)
        real(kind=8)        :: vip(*)
        real(kind=8)        :: sigxp
        real(kind=8)        :: etan
        integer             :: codret
    end subroutine compor_3d_fibre
end interface
