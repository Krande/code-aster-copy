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
interface
    subroutine dpvpre(mod, nvi, option, crit, instam,&
                      instap, nbmat, materf, sigm, deps,&
                      vim, vip, sig, nbre, dsidep,&
                      iret)
        integer(kind=8) :: nbmat
        integer(kind=8) :: nvi
        character(len=8) :: mod
        character(len=16) :: option
        real(kind=8) :: crit(3)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: materf(nbmat, 2)
        real(kind=8) :: sigm(6)
        real(kind=8) :: deps(6)
        real(kind=8) :: vim(nvi)
        real(kind=8) :: vip(nvi)
        real(kind=8) :: sig(6)
        integer(kind=8) :: nbre
        real(kind=8) :: dsidep(6, 6)
        integer(kind=8) :: iret
    end subroutine dpvpre
end interface
