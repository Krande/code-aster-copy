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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nonlinDSConvergenceInit(ds_conv, list_func_acti, ds_contact,&
                                   model)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisr.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/getvr8.h"
#include "asterfort/GetResi.h"
#include "asterfort/SetResi.h"
#include "asterfort/utmess.h"
!
type(NL_DS_Conv), intent(inout) :: ds_conv
integer, optional, intent(in) :: list_func_acti(*)
type(NL_DS_Contact), optional, intent(in) :: ds_contact
character(len=24), optional, intent(in) :: model
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
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    character(len=24) :: sdcont_defi
    character(len=8)  :: exicoq, exipou
    real (kind=8) :: resi_glob_rela, resi_frot, resi_geom,pene_maxi_user
    integer :: iret
    aster_logical :: l_newt_frot, l_newt_geom, l_resi_user, l_rela, l_maxi, l_refe, l_comp
    aster_logical :: l_pena_cont, l_cont, l_geom_sans
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_7')
    endif
!
! - No information from user: RESI_GLOB_RELA with 1E-6
!
    call GetResi(ds_conv, type = 'RESI_GLOB_RELA' , user_para_ = resi_glob_rela,&
                 l_resi_test_ = l_rela)
    call GetResi(ds_conv, type = 'RESI_GLOB_MAXI' , l_resi_test_ = l_maxi)
    call GetResi(ds_conv, type = 'RESI_REFE_RELA' , l_resi_test_ = l_refe)
    call GetResi(ds_conv, type = 'RESI_COMP_RELA' , l_resi_test_ = l_comp)
    l_resi_user = l_rela .or. l_maxi .or. l_refe .or. l_comp
    if (.not.l_resi_user) then
        call SetResi(ds_conv   , type_ = 'RESI_GLOB_RELA', &
                     user_para_ = 1.d-6, l_resi_test_ = ASTER_TRUE)
    endif
!
! - RESI_REFE_RELA with shell and beam
!
    if (l_refe .and. present(model))then
        call dismoi('EXI_COQUE',  model, 'MODELE', repk=exicoq)
        call dismoi('EXI_POUTRE',  model, 'MODELE', repk=exipou)
        if (exicoq .eq. 'OUI' .and. exipou .eq. 'OUI')then
            call utmess('A', 'MECANONLINE5_32')
        elseif (exicoq .eq. 'OUI') then
            call utmess('A', 'MECANONLINE5_38')
        endif
    endif
!
! - Relaxation of convergence criterion: alarm !
!
    if (l_rela .and. resi_glob_rela .gt. 1.0001d-4) then
        call utmess('A', 'MECANONLINE5_21')
    endif
!
! - ARRET=NON: alarm !
!
    if (.not.ds_conv%l_stop) then
        call utmess('A', 'MECANONLINE5_37')
    endif
!
! - No NEWTON/PAS_MINI_ELAS parameter => using ITER_GLOB_MAXI instead of ITER_GLOB_ELAS
!
    call getvr8('NEWTON', 'PAS_MINI_ELAS', iocc=1, nbret=iret)
    if (iret.eq.0) then
        ds_conv%iter_glob_elas = ds_conv%iter_glob_maxi
    endif
!
! - Set contact residuals (not for SIMU_POINT_MAT)
!
    if (present(list_func_acti)) then
!
        l_cont      = isfonc(list_func_acti,'CONTACT')
        sdcont_defi = ds_contact%sdcont_defi
        if (l_cont) then
            l_geom_sans = cfdisl(ds_contact%sdcont_defi, 'REAC_GEOM_SANS')
            if (.not.ds_conv%l_stop .and. .not. l_geom_sans) then
                call utmess('A', 'MECANONLINE5_54')
            endif
        endif
!
! ----- Active functionnalites
!
        l_newt_frot = isfonc(list_func_acti,'FROT_NEWTON')
        l_newt_geom = isfonc(list_func_acti,'GEOM_NEWTON')
        l_pena_cont = isfonc(list_func_acti,'EXIS_PENA')
!
! ----- Activation of contact residuals for generalized Newton
!
        if (l_newt_frot) then
            resi_frot = cfdisr(sdcont_defi, 'RESI_FROT')
            call SetResi(ds_conv   , type_ = 'RESI_FROT', user_para_ = resi_frot,&
                         l_resi_test_ = ASTER_TRUE)
        endif
        if (l_newt_geom) then
            resi_geom = cfdisr(sdcont_defi, 'RESI_GEOM')
            call SetResi(ds_conv   , type_ = 'RESI_GEOM', user_para_ = resi_geom,&
                         l_resi_test_ = ASTER_TRUE)
        endif
        if (l_pena_cont) then
            pene_maxi_user = cfdisr(sdcont_defi, 'PENE_MAXI')
            ! Attention ce parametre est multiplie par la plus petite maille de la zone maitre
            !courante
            ! dans mmalgo
            call SetResi(ds_conv   , type_ = 'RESI_PENE', user_para_ = pene_maxi_user,&
                         l_resi_test_ = ASTER_TRUE)
        endif

    endif
!
! - For line search
!
    ds_conv%line_sear_coef = 1.d0
    ds_conv%line_sear_iter = 0
!
end subroutine
