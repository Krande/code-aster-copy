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
    subroutine vechme(stop     , modelz, lload_namez, lload_infoz, inst        ,&
                      cara_elem, mate  , mateco , vect_elemz , varc_currz , ligrel_calcz,&
                      nharm, basez)
        character(len=1), intent(in) :: stop
        character(len=*), intent(in) :: modelz
        character(len=*), intent(in) :: lload_namez
        character(len=*), intent(in) :: lload_infoz
        real(kind=8), intent(in) :: inst(3)
        character(len=*), intent(in) :: cara_elem
        character(len=*), intent(in) :: mate, mateco
        character(len=*), intent(inout) :: vect_elemz
        character(len=*), optional, intent(in) :: varc_currz
        character(len=*), optional, intent(in) :: ligrel_calcz
        integer, optional, intent(in) :: nharm
        character(len=1), optional, intent(in) :: basez
    end subroutine vechme
end interface
