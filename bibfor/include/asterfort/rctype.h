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
    subroutine rctype(jmat     , nb_para_list, para_list_name, para_list_vale, para_vale,&
                      para_type, keyw_factz  , keywz, materi)
        integer(kind=8), intent(in) :: jmat
        integer(kind=8), intent(in) :: nb_para_list
        character(len=*), intent(in) :: para_list_name(*)
        real(kind=8), intent(in) :: para_list_vale(*)
        real(kind=8), intent(out) :: para_vale
        character(len=*), intent(out) :: para_type
        character(len=*), optional, intent(in) :: keyw_factz
        character(len=*), optional, intent(in) :: keywz
        character(len=*), optional, intent(in) :: materi
    end subroutine rctype
end interface
