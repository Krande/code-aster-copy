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
!
subroutine nonlinDSAlgoParaInit(list_func_acti, ds_algopara, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/cfdisl.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm parameters management
!
! Initializations for algorithm parameters management
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! IO  ds_algopara      : datastructure for algorithm parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: reli_rho_mini, reli_rho_maxi, reli_rho_excl, swap
    aster_logical :: l_cont_disc, l_unil, l_unil_pena
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_6')
    end if
!
! - Active functionnalites
!
    l_cont_disc = isfonc(list_func_acti, 'CONT_DISCRET')
    l_unil = isfonc(list_func_acti, 'LIAISON_UNILATER')
!
! - Symmetric matrix for DISCRETE contact
!
    if (l_cont_disc) then
        ds_algopara%l_matr_rigi_syme = .true._1
    end if
    if (l_unil) then
        l_unil_pena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
        if (.not. l_unil_pena) then
            ds_algopara%l_matr_rigi_syme = .true._1
        end if
    end if
!
! - Update bounds for line search
!
    if (ds_algopara%l_line_search) then
        reli_rho_mini = ds_algopara%line_search%rho_mini
        reli_rho_maxi = ds_algopara%line_search%rho_maxi
        reli_rho_excl = ds_algopara%line_search%rho_excl
        if (reli_rho_mini .ge. -reli_rho_excl .and. reli_rho_mini .le. reli_rho_excl) then
            call utmess('A', 'MECANONLINE5_46')
            reli_rho_mini = +reli_rho_excl
        end if
        if (reli_rho_maxi .ge. -reli_rho_excl .and. reli_rho_maxi .le. reli_rho_excl) then
            call utmess('A', 'MECANONLINE5_47')
            reli_rho_maxi = -reli_rho_excl
        end if
        if (reli_rho_maxi .lt. reli_rho_mini) then
            call utmess('A', 'MECANONLINE5_44')
            swap = reli_rho_maxi
            reli_rho_maxi = reli_rho_mini
            reli_rho_mini = swap
        end if
        if (abs(reli_rho_maxi-reli_rho_mini) .le. r8prem()) then
            call utmess('F', 'MECANONLINE5_43')
        end if
        ds_algopara%line_search%rho_mini = reli_rho_mini
        ds_algopara%line_search%rho_maxi = reli_rho_maxi
    end if
!
end subroutine
