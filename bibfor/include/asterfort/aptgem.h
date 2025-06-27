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
    subroutine aptgem(sdappa , mesh     , newgeo   , sdcont_defi, model_ndim,&
                      i_zone , zone_type, epsi_maxi  , jdecma    ,&
                      nb_elem, err_appa)
        character(len=19), intent(in) :: sdappa
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: sdcont_defi
        character(len=19), intent(in) :: newgeo
        integer(kind=8), intent(in) :: model_ndim
        integer(kind=8), intent(in) :: i_zone
        integer(kind=8), intent(in) :: jdecma
        integer(kind=8), intent(in) :: nb_elem
        character(len=4), intent(in) :: zone_type
        real(kind=8), intent(in) :: epsi_maxi
        integer(kind=8), intent(inout) :: err_appa
    end subroutine aptgem
end interface
