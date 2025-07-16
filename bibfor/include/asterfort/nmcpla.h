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
    subroutine nmcpla(BEHinteg,&
                      fami, kpg, ksp, ndim, typmod, imat, &
                      compor_plas, compor_creep, carcri, &
                      timed, timef, neps, epsdt, depst, &
                      nsig, sigd, vind, option, &
                      sigf, vinf, ndsde, dsde, iret)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer(kind=8) :: ndsde
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        integer(kind=8) :: ndim
        character(len=8) :: typmod(*)
        integer(kind=8) :: imat
        character(len=16), intent(in) :: compor_plas(*)
        character(len=16), intent(in) :: compor_creep(*)
        real(kind=8), intent(in) :: carcri(*)
        real(kind=8) :: timed
        real(kind=8) :: timef
        integer(kind=8) :: neps
        real(kind=8) :: epsdt(6)
        real(kind=8) :: depst(6)
        integer(kind=8) :: nsig
        real(kind=8) :: sigd(6)
        real(kind=8) :: vind(*)
        character(len=16) :: option
        real(kind=8) :: sigf(6)
        real(kind=8) :: vinf(*)
        real(kind=8) :: dsde(ndsde)
        integer(kind=8) :: iret
    end subroutine nmcpla
end interface
