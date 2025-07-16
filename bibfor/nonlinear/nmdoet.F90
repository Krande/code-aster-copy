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
subroutine nmdoet(model, compor, list_func_acti, nume_ddl, sdpilo, &
                  sddyna, ds_errorindic, hval_algo, l_acce_zero, ds_inout, &
                  ds_energy)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndloam.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdoin.h"
#include "asterfort/nmetl1.h"
#include "asterfort/nmetl2.h"
#include "asterfort/nmetl3.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/infdbg.h"
#include "asterfort/nonlinDSEnergyInitValues.h"
!
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
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Init
!
! Read initial state
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  compor           : name of <CARTE> COMPOR
! In  list_func_acti   : list of active functionnalities
! In  nume_ddl         : name of nume_ddl object (numbering equation)
! In  sddyna           : dynamic parameters datastructure
! In  hval_algo        : hat-variable for algorithms fields
! In  sdpilo           : continuation ("PILOTAGE") parameters datastructure
! IO  ds_errorindic    : datastructure for error indicator
! IO  ds_inout         : datastructure for input/output management
! Out l_acce_zero      : .true. if initial acceleration must been computed
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_field
    character(len=24) :: field_type
    aster_logical :: l_stin_evol, l_state_init
    integer(kind=8) :: nb_equa, init_nume, iret, i, i_field
    character(len=8) :: calcri, stin_evol
    character(len=24) :: typpil, typsel
    character(len=19) :: depold
    character(len=24) :: champ1, champ2, dep2, dep1
    integer(kind=8) :: jv_para
    aster_logical :: l_pilo, lpiarc, l_cont_cont
    aster_logical :: l_expl_gene, l_reuse, l_erre_thm, l_mstp
    aster_logical :: l_zero, l_acti, l_ener, l_read, verbose
    real(kind=8) :: coefav, init_time
    real(kind=8), pointer :: plir(:) => null()
    character(len=24), pointer :: pltk(:) => null()
    real(kind=8), pointer :: vdep1(:) => null()
    real(kind=8), pointer :: vdep2(:) => null()
    real(kind=8), pointer :: depol(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_26')
    end if
    verbose = niv .ge. 2
!
! - Initializations
!
    dep1 = '&&CNPART.CHP1'
    dep2 = '&&CNPART.CHP2'
    l_acce_zero = .false.
    lpiarc = .false.
    call dismoi('NB_EQUA', nume_ddl, 'NUME_DDL', repi=nb_equa)
!
! - Model can compute rigidity ?
!
    call dismoi('CALC_RIGI', model, 'MODELE', repk=calcri)
    if (calcri .ne. 'OUI') then
        call utmess('F', 'CALCULEL2_65', sk=model)
    end if
!
! - Parameters from input/output datastructure
!
    nb_field = ds_inout%nb_field
!
! - Active functionnalities
!
    l_pilo = isfonc(list_func_acti, 'PILOTAGE')
    l_cont_cont = isfonc(list_func_acti, 'CONT_CONTINU')
    l_expl_gene = ndynlo(sddyna, 'EXPL_GENE')
    l_reuse = isfonc(list_func_acti, 'REUSE')
    l_erre_thm = isfonc(list_func_acti, 'ERRE_TEMPS_THM')
    l_ener = isfonc(list_func_acti, 'ENERGIE')
    l_mstp = ndynlo(sddyna, 'MULTI_PAS')
!
! - Get previous displacement field
!
    call nmchex(hval_algo, 'SOLALG', 'DEPOLD', depold)
!
! - Does ETAT_INIT (initial state) exist ?
!
    l_state_init = isfonc(list_func_acti, 'ETAT_INIT')
!
! - PILOTAGE LONGUEUR D'ARC AVEC ANGL_INCR_DEPL: IL FAUT LES DEUX
! - DERNIERS DEPLACEMENTS POUR QUE CA MARCHE (CHAMP DEPOLD)
!
    if (l_pilo) then
        call jeveuo(sdpilo(1:19)//'.PLTK', 'L', vk24=pltk)
        typpil = pltk(1)
        typsel = pltk(6)
        lpiarc = .false.
        if (typpil .eq. 'LONG_ARC') then
            if (typsel(1:14) .eq. 'ANGL_INCR_DEPL') then
                lpiarc = .true.
            end if
        end if
    end if
!
! - Get name of result datastructure in ETAT_INIT
!
    l_stin_evol = ds_inout%l_stin_evol
    stin_evol = ds_inout%stin_evol
!
! - ALARME SI CONTACT CONTINU AVEC UN CONCEPT REENTRANT
!
    if (l_cont_cont) then
        if (l_reuse) then
            if (.not. isfonc(list_func_acti, 'CONTACT_INIT')) then
                call utmess('A', 'MECANONLINE4_14')
            end if
        else if (l_stin_evol) then
            if (.not. isfonc(list_func_acti, 'CONTACT_INIT')) then
                call utmess('A', 'MECANONLINE4_15')
            end if
        end if
    end if
!
! - Initial storing index and time
!
    call nmdoin(ds_inout)
    init_time = ds_inout%init_time
    init_nume = ds_inout%init_nume
!
! - Print
!
    if (l_state_init) then
        call utmess('I', 'ETATINIT_10')
        if (l_stin_evol) then
            call utmess('I', 'ETATINIT_11', sk=stin_evol, sr=init_time, si=init_nume)
        else
            if (l_ener) then
                call utmess('I', 'ETATINIT_5')
            end if
            if (init_nume .eq. -1) then
                call utmess('I', 'ETATINIT_20')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    else
        if (l_reuse) then
            call utmess('A', 'ETATINIT_1')
        else
            call utmess('I', 'ETATINIT_20')
        end if
    end if
!
! - Loop on fields
!
    do i_field = 1, nb_field
!
! ----- Read field for ETAT_INIT - From results datastructure
!
        if (l_stin_evol) then
            call nmetl1(i_field, ds_inout)
        end if
!
! ----- Read field for ETAT_INIT - Field by field
!
        call nmetl2(model, i_field, ds_inout)
!
! ----- Read field for ETAT_INIT - Some checks
!
        call nmetl3(model, compor, i_field, ds_inout, verbose)
    end do
!
! - VERIFICATION COMPATIBILITE PILOTAGE
!
    if (l_stin_evol .and. lpiarc) then
        call rsexch(' ', stin_evol, 'DEPL', init_nume, champ1, iret)
        call rsexch(' ', stin_evol, 'DEPL', init_nume-1, champ2, iret)
        if (iret .ne. 0) then
            call utmess('F', 'MECANONLINE4_47', sk=stin_evol)
        end if
        call vtcopy(champ1, dep1, iret)
        if (iret .ne. 0) then
            call utmess("F", "FIELD0_9")
        end if
        call vtcopy(champ2, dep2, iret)
        if (iret .ne. 0) then
            call utmess("F", "FIELD0_9")
        end if
        call jeveuo(dep1(1:19)//'.VALE', 'L', vr=vdep1)
        call jeveuo(dep2(1:19)//'.VALE', 'L', vr=vdep2)
        call jeveuo(depold(1:19)//'.VALE', 'E', vr=depol)
        do i = 1, nb_equa
            depol(i) = vdep1(i)-vdep2(i)
        end do
        call jeveuo(sdpilo(1:19)//'.PLIR', 'E', vr=plir)
        call rsadpa(stin_evol, 'L', 1, 'COEF_MULT', init_nume, &
                    0, sjv=jv_para, istop=0)
        coefav = zr(jv_para)
        if (coefav .ne. 0.d0 .and. coefav .ne. r8vide()) then
            plir(6) = coefav
        end if
    end if
!
! - LECTURE DES INDICATEURS D'ERREUR EN TEMPS EN THM
!
    if (l_stin_evol .and. l_erre_thm) then
        call rsadpa(stin_evol, 'L', 1, 'ERRE_TPS_LOC', init_nume, &
                    0, sjv=jv_para, istop=0)
        ds_errorindic%erre_thm_loca = zr(jv_para)
        call rsadpa(stin_evol, 'L', 1, 'ERRE_TPS_GLOB', init_nume, &
                    0, sjv=jv_para, istop=0)
        ds_errorindic%erre_thm_glob = zr(jv_para)
    end if
!
! - lecture des données initiales sur l'énergie
!
    if (l_stin_evol .and. l_ener) then
        call nonlinDSEnergyInitValues(ds_energy, stin_evol, ds_inout)
    end if
!
! - CAS DE LA DYNAMIQUE: VITESSE ET ACCELERATION INITIALES
!
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        l_zero = ds_inout%field(i_field)%init_type .eq. 'ZERO'
        l_acti = ds_inout%l_field_acti(i_field)
        if (field_type .eq. 'VITE') then
            if (l_zero .and. l_acti) then
                call utmess('I', 'MECANONLINE4_22')
            end if
        end if
        if (field_type .eq. 'ACCE') then
            if (l_zero .and. l_acti) then
                call utmess('I', 'MECANONLINE4_23')
                l_acce_zero = .true.
            end if
        end if
    end do

!
!   Check initial stress state for multi-step schemes
!
    if (l_mstp) then
        do i_field = 1, nb_field
            field_type = ds_inout%field(i_field)%type
            l_acti = ds_inout%l_field_acti(i_field)
            l_read = ds_inout%field(i_field)%init_type .eq. 'READ'
            if ((field_type .eq. 'SIEF_ELGA') .and. (.not. l_stin_evol)) then
                if (l_acti .and. l_read) then
                    call utmess('A', 'MECANONLINE4_50')
                end if
            end if
        end do
    end if

!
! - PROJECTION MODALE EN EXPLICITE
!
    if (l_expl_gene) then
        call ndloam(sddyna, stin_evol, l_stin_evol, init_nume)
    end if
!
    call jedema()
end subroutine
