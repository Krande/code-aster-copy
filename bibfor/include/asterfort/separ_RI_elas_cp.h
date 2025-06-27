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
    subroutine separ_RI_elas_cp(elas_id ,nu , g, nui ,gi, &
                             e1     , e2  ,&
                             nu12   , &
                             e1i     , e2i  ,&
                             nu12i   , &
                             hr, hi)
!
        integer(kind=8), intent(in) :: elas_id
        real(kind=8), intent(in) :: nu, g, e1, e2
        real(kind=8), intent(in) :: nu12
        real(kind=8), intent(in) :: nui, gi, e1i, e2i
        real(kind=8), intent(in) :: nu12i
        real(kind=8), intent(out) :: hr(3), hi(3)
    end subroutine separ_RI_elas_cp
end interface
