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
    subroutine lcsena(elem_dime, nb_lagr, nb_node_slav, indi_lagc, &
                      lagrc    , vtmp)
        integer(kind=8), intent(in) :: elem_dime
        integer(kind=8), intent(in) :: nb_lagr
        integer(kind=8), intent(in) :: nb_node_slav
        integer(kind=8), intent(in) :: indi_lagc(10)
        real(kind=8), intent(in) :: lagrc
        real(kind=8), intent(inout) :: vtmp(55)
    end subroutine lcsena
end interface
