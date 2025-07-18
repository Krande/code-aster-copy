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
#include "asterf_types.h"
!
interface
    subroutine nmdoet(model , compor       , list_func_acti, nume_ddl   , sdpilo  ,&
                      sddyna, ds_errorindic, hval_algo     , l_acce_zero, ds_inout,&
                      ds_energy)
        use NonLin_Datastructure_type
        character(len=24), intent(in) :: model
        character(len=24), intent(in) :: compor
        type(NL_DS_ErrorIndic), intent(inout) :: ds_errorindic
        character(len=24), intent(in) :: nume_ddl
        character(len=19), intent(in) :: sddyna
        character(len=19), intent(in) :: sdpilo
        character(len=19), intent(in) :: hval_algo(*)
        integer(kind=8), intent(in) :: list_func_acti(*)
        aster_logical, intent(out) :: l_acce_zero
        type(NL_DS_InOut), intent(inout) :: ds_inout
        type(NL_DS_Energy), intent(inout) :: ds_energy
    end subroutine nmdoet
end interface
