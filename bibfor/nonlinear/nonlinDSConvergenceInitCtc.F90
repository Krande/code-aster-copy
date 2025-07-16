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
subroutine nonlinDSConvergenceInitCtc(ds_conv, list_func_acti, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisr.h"
#include "asterfort/cfdisl.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/SetResi.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_Conv), intent(inout) :: ds_conv
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Initializations for convergence management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_conv          : datastructure for convergence management
! In  list_func_acti   : list of active functionnalities
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=24) :: sdcont_defi
    real(kind=8) ::resi_frot, resi_geom, pene_maxi_user
    aster_logical :: l_newt_frot, l_newt_geom
    aster_logical :: l_pena_cont, l_cont, l_geom_sans
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)

! - Active functionnalites
    l_cont = isfonc(list_func_acti, 'CONTACT')
    l_newt_frot = isfonc(list_func_acti, 'FROT_NEWTON')
    l_newt_geom = isfonc(list_func_acti, 'GEOM_NEWTON')
    l_pena_cont = isfonc(list_func_acti, 'EXIS_PENA')

! - Access to datastructture
    sdcont_defi = ds_contact%sdcont_defi
    if (l_cont) then
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_13')
        end if
    end if

! - Alarm when using ARRET='NON' with contact
    if (l_cont) then
        l_geom_sans = cfdisl(ds_contact%sdcont_defi, 'REAC_GEOM_SANS')
        if (.not. ds_conv%l_stop .and. .not. l_geom_sans) then
            call utmess('A', 'MECANONLINE5_54')
        end if
    end if

! - Activation of contact residuals for generalized Newton
    if (l_newt_frot) then
        resi_frot = cfdisr(sdcont_defi, 'RESI_FROT')
        call SetResi(ds_conv, type_='RESI_FROT', user_para_=resi_frot, &
                     l_resi_test_=ASTER_TRUE)
    end if
    if (l_newt_geom) then
        resi_geom = cfdisr(sdcont_defi, 'RESI_GEOM')
        call SetResi(ds_conv, type_='RESI_GEOM', user_para_=resi_geom, &
                     l_resi_test_=ASTER_TRUE)
    end if
    if (l_pena_cont) then
        pene_maxi_user = cfdisr(sdcont_defi, 'PENE_MAXI')
! Attention ce parametre est multiplie par la plus petite maille de la zone maitre
! courante
! dans mmalgo
        call SetResi(ds_conv, type_='RESI_PENE', user_para_=pene_maxi_user, &
                     l_resi_test_=ASTER_TRUE)
    end if
!
end subroutine
