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
#include "asterf_types.h"
!
interface
    subroutine nmeteo(result, sddisc , ds_inout , force, nume_store,&
                      time  , i_field, ds_print_)
        use NonLin_Datastructure_type
        type(NL_DS_InOut), intent(in) :: ds_inout
        character(len=19), intent(in) :: sddisc
        character(len=8), intent(in) :: result
        integer(kind=8), intent(in) :: i_field
        integer(kind=8), intent(in) :: nume_store
        real(kind=8), intent(in) :: time
        aster_logical, intent(in) :: force
        type(NL_DS_Print), optional, intent(in) :: ds_print_
    end subroutine nmeteo
end interface
