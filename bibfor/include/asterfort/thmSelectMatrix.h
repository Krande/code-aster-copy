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
    subroutine thmSelectMatrix(ds_thm,&
                               ndim  , dimdef, inte_type,&
                               addeme, addete, addep1   , addep2, adde2nd, &
                               a     , as    ,&
                               c     , cs    )
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8), intent(in) :: ndim, dimdef
        character(len=3), intent(in) :: inte_type
        integer(kind=8), intent(in) :: addeme, addete, addep1, addep2, adde2nd
        real(kind=8), intent(out) :: a(2), as(2)
        real(kind=8), intent(out) :: c(21), cs(21)
    end subroutine thmSelectMatrix
end interface
