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

subroutine nmdoch_wrap(list_load0, l_load_user0, l_calc_user0, list_load_resu0, base, &
                       ligrel_slav0, ligrel_cont0)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmdoch.h"
!
! person_in_charge: nicolas.sellenet at edf.fr
!
    character(len=*), intent(in) :: list_load0
    integer, intent(in) :: l_load_user0, l_calc_user0
    character(len=*), intent(in) :: list_load_resu0, ligrel_slav0, ligrel_cont0
    character(len=*), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics - Read parameters
!
! Get loads information and create datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  list_load_resu : name of datastructure for list of loads from result datastructure
! In  l_load_user    : .true. if loads come from user (EXCIT)
! In  list_load      : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: list_load
    aster_logical :: l_load_user, l_calc_user
    character(len=19) :: list_load_resu, ligrel_slav, ligrel_cont
!
! --------------------------------------------------------------------------------------------------
!
    l_load_user = int_to_logical(l_load_user0)
    l_calc_user = int_to_logical(l_calc_user0)
    list_load = list_load0
    list_load_resu = list_load_resu0
    ligrel_slav = ligrel_slav0
    ligrel_cont = ligrel_cont0
    call nmdoch(list_load, l_load_user, list_load_resu, base, l_calc_user, &
                ligrel_slav, ligrel_cont)
end subroutine
