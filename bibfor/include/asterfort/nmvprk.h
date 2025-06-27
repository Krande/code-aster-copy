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
interface
    subroutine nmvprk(fami, kpg, ksp, ndim, typmod,&
                      imat, comp, crit, timed, timef,&
                      neps, epsdt, depst, sigd, nvi, vind,&
                      opt, angmas, sigf, vinf, dsde,&
                      iret, mult_comp_)
        integer(kind=8) :: neps
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: ndim
        character(len=8) :: typmod(*)
        integer(kind=8) :: imat
        character(len=16) :: comp(*)
        real(kind=8) :: crit(*)
        real(kind=8) :: timed
        real(kind=8) :: timef
        real(kind=8) :: epsdt(neps)
        real(kind=8) :: depst(neps)
        real(kind=8) :: sigd(6)
        integer(kind=8), intent(in):: nvi
        real(kind=8) :: vind(*)
        character(len=16) :: opt
        real(kind=8) :: angmas(*)
        real(kind=8) :: sigf(6)
        real(kind=8) :: vinf(*)
        real(kind=8) :: dsde(6, *)
        integer(kind=8) :: iret
        character(len=16), optional, intent(in) :: mult_comp_
    end subroutine nmvprk
end interface
