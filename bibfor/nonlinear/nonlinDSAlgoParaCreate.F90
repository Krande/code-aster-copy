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
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1403
!
subroutine nonlinDSAlgoParaCreate(ds_algopara)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(NL_DS_AlgoPara), intent(out) :: ds_algopara
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm parameters management
!
! Create algorithm parameters datastructure
!
! --------------------------------------------------------------------------------------------------
!
! Out ds_algopara      : datastructure for algorithm parameters
!
! --------------------------------------------------------------------------------------------------
!
    ds_algopara%method = ' '
    ds_algopara%matrix_pred = ' '
    ds_algopara%matrix_corr = ' '
    ds_algopara%reac_incr = 0
    ds_algopara%reac_iter = 0
    ds_algopara%pas_mini_elas = -9999.0d0
    ds_algopara%reac_iter_elas = 0
    ds_algopara%l_dyna = ASTER_FALSE
    ds_algopara%l_line_search = ASTER_FALSE
    ds_algopara%l_pilotage = ASTER_FALSE
    ds_algopara%result_prev_disp = ' '
    ds_algopara%l_matr_rigi_syme = ASTER_FALSE
!
! - Parameters for line search
!
    ds_algopara%line_search%method = ' '
    ds_algopara%line_search%resi_rela = 0.d0
    ds_algopara%line_search%iter_maxi = 0
    ds_algopara%line_search%rho_mini = 0.d0
    ds_algopara%line_search%rho_maxi = 0.d0
    ds_algopara%line_search%rho_excl = 0.d0
!
end subroutine
