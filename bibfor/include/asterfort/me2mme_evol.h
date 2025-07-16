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
    subroutine me2mme_evol(modelZ, caraElemZ, mateZ, matecoZ, nharm, jvBase, &
                           iLoad, loadName, ligrel_calcZ, inst_prev, inst_curr, &
                           inst_theta, resuElem, vectElem)
        character(len=*), intent(in) :: modelZ
        character(len=*), intent(in) :: caraElemZ
        character(len=*), intent(in) :: mateZ
        character(len=*), intent(in) :: matecoZ
        integer(kind=8), intent(in) :: nharm
        character(len=1), intent(in) :: jvBase
        integer(kind=8), intent(in) :: iLoad
        character(len=8), intent(in) :: loadName
        character(len=*), intent(in) :: ligrel_calcZ
        real(kind=8), intent(in) :: inst_prev
        real(kind=8), intent(in) :: inst_curr
        real(kind=8), intent(in) :: inst_theta
        character(len=19), intent(inout) :: resuElem
        character(len=19), intent(in) :: vectElem
    end subroutine me2mme_evol
end interface
