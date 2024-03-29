! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine licrsd_save(list_load)
!
    implicit none
!
#include "asterfort/gnomsd.h"
!
    character(len=19), intent(out) :: list_load
!
! --------------------------------------------------------------------------------------------------
!
! List of loads - Utility
!
! Generate name of list of loads to save in results datastructure
!
! --------------------------------------------------------------------------------------------------
!
! Out list_load_out     : list of loads to save
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: noobj
!
! --------------------------------------------------------------------------------------------------
!
    list_load = ' '
    noobj = '12345678'//'.1234'//'.EXCIT'
    call gnomsd(' ', noobj, 10, 13)
    list_load = noobj(1:19)
!
end subroutine
