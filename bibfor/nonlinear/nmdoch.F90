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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmdoch(list_load, l_load_user, list_load_resu_, base, l_calc_user, &
                  ligrel_slav, ligrel_cont)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmdoch_nbload.h"
#include "asterfort/load_list_getp.h"
#include "asterfort/load_unde_diri.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/focste.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liscad.h"
#include "asterfort/lisccr.h"
#include "asterfort/lisexp.h"
#include "asterfort/lislfc.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/loadGetNeumannType.h"
#include "asterfort/loadExcludedForAnalysis.h"
!
    character(len=19), intent(in) :: list_load
    aster_logical, intent(in) :: l_load_user
    aster_logical, optional, intent(in) :: l_calc_user
    character(len=19), optional, intent(in) :: list_load_resu_, ligrel_slav, ligrel_cont
    character(len=1), optional, intent(in) :: base
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
! In  l_load_user    : .true. if loads come from user (command file)
! In  list_load      : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nb_info_maxi = 99
    character(len=24) :: list_info_type(nb_info_maxi)
    integer :: n1, npilo, nb_load, nb_load_cont
    character(len=1) :: bas
    integer ::  i_excit, i_load, iret, i_load_new
    real(kind=8) :: coef
    character(len=4) :: typcal
    character(len=8) :: k8bid, load_type, func_para_inst, const_func
    character(len=16) :: nomcmd, typesd, load_apply, load_keyword
    character(len=8) :: load_name, load_func, model, load_model
    character(len=24) :: info_type
    character(len=19) :: lisdbl, list_load_resu, func_cont
    character(len=24) :: ligrch, lchin
    integer :: i_diri_suiv
    aster_logical :: l_func_c, l_zero_allowed, l_diri_undead, l_stat, l_calc
    integer :: nb_info_type
    integer, pointer :: v_ll_infc(:) => null()
    integer, pointer :: v_llresu_info(:) => null()
    character(len=24), pointer :: v_llresu_name(:) => null()
    character(len=8), pointer :: v_list_dble(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call getres(k8bid, typesd, nomcmd)
    bas = 'V'
    if (present(base)) then
        bas = base
    end if
!
! - Initializations
!
    l_stat = nomcmd .eq. 'STAT_NON_LINE'
    list_load_resu = ' '
    if (present(list_load_resu_)) then
        list_load_resu = list_load_resu_
    end if
    nb_load = 0
    nb_load_cont = 0
    const_func = '&&NMDOME'
    lisdbl = '&&NMDOME.LISDBL'
    l_func_c = ASTER_FALSE
    info_type = 'RIEN'
    npilo = 0
    i_excit = 0
    i_load_new = 0
    i_diri_suiv = 0
    l_diri_undead = ASTER_FALSE
!
! - Command for computation ?
!
    if (present(l_calc_user)) then
        l_calc = l_calc_user
    else
        l_calc = nomcmd .eq. 'DYNA_NON_LINE' .or. &
                 nomcmd .eq. 'STAT_NON_LINE' .or. &
                 nomcmd .eq. 'MECA_STATIQUE' .or. &
                 nomcmd .eq. 'CALC_FORC_NONL' .or. &
                 nomcmd .eq. 'DYNA_VIBRA'
    end if
!
! - Get model
!
    model = ' '
    if (nomcmd .ne. 'CREA_RESU') then
        call getvid(' ', 'MODELE', scal=model, nbret=n1)
    end if
!
! - Can we create "zero-load" list of loads datastructure ?
!
    l_zero_allowed = ASTER_FALSE
    if (l_load_user) then
        l_zero_allowed = nomcmd .eq. 'DYNA_NON_LINE' .or. &
                         nomcmd .eq. 'STAT_NON_LINE' .or. &
                         nomcmd .eq. 'MECA_STATIQUE' .or. &
                         nomcmd .eq. 'CALC_FORC_NONL' .or. &
                         nomcmd .eq. 'CALC_CHAMP' .or. &
                         nomcmd .eq. 'CALCUL' .or. &
                         nomcmd .eq. 'POST_ELEM' .or. &
                         nomcmd .eq. 'LIRE_RESU' .or. &
                         nomcmd .eq. 'CREA_RESU'
    else
        l_zero_allowed = ASTER_TRUE
    end if
!
! - Get number of loads for loads datastructure
!
    call nmdoch_nbload(l_load_user, list_load_resu, l_zero_allowed, nb_load, &
                       load_keyword)
!
! - Get number of contact ligrel for loads datastructure
!
    if (present(ligrel_slav) .and. ligrel_slav .ne. ' ') then
        nb_load_cont = nb_load_cont+1
    end if
    if (present(ligrel_cont) .and. ligrel_cont .ne. ' ') then
        nb_load_cont = nb_load_cont+1
    end if
!
! - Create "zero-load" list of loads datastructure
!
    if (nb_load+nb_load_cont .eq. 0) then
        call lisccr('MECA', list_load, 1, bas)
        call jeveuo(list_load(1:19)//'.INFC', 'E', vi=v_ll_infc)
        v_ll_infc(1) = 1
    end if
!
! - Access to saved list of loads datastructure
!
    if (.not. l_load_user) then
        call jeveuo(list_load_resu(1:19)//'.INFC', 'L', vi=v_llresu_info)
        call jeveuo(list_load_resu(1:19)//'.LCHA', 'L', vk24=v_llresu_name)
    end if
!
    if (nb_load+nb_load_cont .ne. 0) then
        ASSERT(load_keyword .ne. 'None')
!
! ----- Create list of loads
!
        call lisccr('MECA', list_load, nb_load+nb_load_cont, bas)
!
! ----- List of loads to avoid same loads
!
        AS_ALLOCATE(vk8=v_list_dble, size=nb_load)
!
! ----- Loop on loads
!
        do i_load = 1, nb_load
!
            if (.not. l_load_user) then
                if (v_llresu_name(i_load) .eq. ' ') cycle
            end if
!
! --------- Get parameters for construct list of loads
!
            call load_list_getp('MECA', l_load_user, v_llresu_info, v_llresu_name, &
                                v_list_dble, load_keyword, i_load, nb_load, &
                                i_excit, load_name, load_type, ligrch, &
                                load_apply)
!
! --------- Check model / model of load
!
            if (n1 .eq. 1 .and. model .ne. ' ') then
                call dismoi('NOM_MODELE', load_name, 'CHARGE', repk=load_model)
                if (load_model .ne. model) then
                    call utmess('F', 'CHARGES_33', nk=3, valk=[load_name, load_model, model])
                end if
            end if
!
! --------- Get function applied to load
!
            if (nomcmd(1:10) .eq. 'DYNA_VIBRA') then
                call getvtx(' ', 'TYPE_CALCUL', iocc=1, scal=typcal)
                l_func_c = (typcal .eq. 'HARM')
            else
                l_func_c = (nomcmd .eq. 'LIRE_RESU' .and. typesd .eq. 'DYNA_HARMO')
            end if
            call lislfc(list_load_resu, i_load, i_excit, l_load_user, &
                        l_func_c, load_keyword, const_func, load_func, bas)
!
! --------- Check loads "PILOTAGE"
!
            if (load_apply .eq. 'FIXE_PILO' .or. load_apply .eq. 'SUIV_PILO') then
                npilo = npilo+1
                if (load_func .ne. const_func) then
                    call utmess('F', 'CHARGES_38', sk=load_name)
                end if
            end if
!
! --------- Dirichlet loads (AFFE_CHAR_CINE)
!
            nb_info_type = 0
            info_type = 'RIEN'
            if (load_type(1:5) .eq. 'CIME_') then
                if (load_apply .eq. 'SUIV') then
                    call utmess('F', 'CHARGES_23', sk=load_name)
                else if (load_apply .eq. 'FIXE_PILO') then
                    call utmess('F', 'CHARGES_27', sk=load_name)
                else if (load_apply .eq. 'DIDI') then
                    call utmess('F', 'CHARGES_24', sk=load_name)
                else if (load_apply .eq. 'FIXE_CSTE') then
                    if (load_type(5:7) .eq. '_FT') then
                        info_type = 'CINE_FT'
                    else if (load_type(5:7) .eq. '_FO') then
                        info_type = 'CINE_FO'
                    else
                        info_type = 'CINE_CSTE'
                    end if
                else
                    ASSERT(.false.)
                end if
            end if
            if (info_type .ne. 'RIEN') then
                nb_info_type = nb_info_type+1
                ASSERT(nb_info_type .lt. nb_info_maxi)
                list_info_type(nb_info_type) = info_type
            end if
!
! --------- Dirichlet loads (AFFE_CHAR_MECA)
!
            info_type = 'RIEN'
            lchin = ligrch(1:13)//'.CIMPO.DESC'
            call jeexin(lchin, iret)
            if (iret .ne. 0) then
                call dismoi('PARA_INST', lchin, 'CARTE', repk=func_para_inst)
                if (load_apply .eq. 'SUIV') then
                    info_type = 'DIRI_SUIV'
                    l_diri_undead = .true.
                else if (load_apply .eq. 'FIXE_PILO') then
                    if (load_type(5:7) .eq. '_FO') then
                        info_type = 'DIRI_PILO_F'
                        if (func_para_inst(1:3) .eq. 'OUI') then
                            call utmess('F', 'CHARGES_28', sk=load_name)
                        end if
                    else
                        info_type = 'DIRI_PILO'
                    end if
                else if (load_apply .eq. 'DIDI') then
                    if (load_type(5:7) .eq. '_FO') then
                        info_type = 'DIRI_FO_DIDI'
                        if (func_para_inst(1:3) .eq. 'OUI') then
                            info_type = 'DIRI_FT_DIDI'
                        end if
                    else
                        info_type = 'DIRI_CSTE_DIDI'
                    end if
                else if (load_apply .eq. 'FIXE_CSTE') then
                    if (load_type(5:7) .eq. '_FO') then
                        info_type = 'DIRI_FO'
                        if (func_para_inst(1:3) .eq. 'OUI') then
                            info_type = 'DIRI_FT'
                        end if
                    else
                        info_type = 'DIRI_CSTE'
                    end if
                else
                    ASSERT(.false.)
                end if
            end if
            if (info_type .ne. 'RIEN') then
                nb_info_type = nb_info_type+1
                ASSERT(nb_info_type .lt. nb_info_maxi)
                list_info_type(nb_info_type) = info_type
            end if
!
! --------- Get Neuman loads
!
            call loadGetNeumannType(l_stat, load_name, ligrch, &
                                    load_apply, load_type, &
                                    nb_info_type, nb_info_maxi, list_info_type)
!
! --------- Add new load(s) in list
!

            if (nb_info_type .gt. 0) then
                i_load_new = i_load_new+1
                call liscad('MECA', list_load, i_load_new, load_name, load_func, &
                            nb_info_type, list_info_type)
            end if
!
        end do
!
! ---- Add contact ligrel
!
        if (nb_load_cont > 0) then
!
! --- Prepare constant function
!
            func_cont = '&&NMDOME.FCONT'
            call jeexin(func_cont//'.PROL', iret)
            if (iret .eq. 0) then
                coef = 1.d0
                call focste(func_cont, 'TOUTRESU', coef, bas)
            end if

            if (present(ligrel_slav) .and. ligrel_slav .ne. ' ') then
                i_load_new = i_load_new+1
                call liscad('MECA', list_load, i_load_new, ligrel_slav, func_cont, &
                            info_typez='ELEM_TARDIF')
            end if
            if (present(ligrel_cont) .and. ligrel_cont .ne. ' ') then
                i_load_new = i_load_new+1
                call liscad('MECA', list_load, i_load_new, ligrel_cont, func_cont, &
                            info_typez='ELEM_TARDIF')
            end if
        end if
!
! ---- PILOTAGE POSSIBLE SI IL YA DES CHARGES PILOTEES !
!
        if (nomcmd .ne. 'LIRE_RESU') then
            if (nomcmd .eq. 'STAT_NON_LINE') then
                call getvtx('PILOTAGE', 'TYPE', iocc=1, nbret=n1)
                if (n1 .ne. 0 .and. npilo .eq. 0) then
                    call utmess('F', 'CHARGES_39')
                end if
                if (npilo .gt. 1) then
                    call utmess('F', 'CHARGES_40')
                end if
            end if
        end if
!
! ----- Some loads are prohibited with PILOTAGE
!
        if (npilo .ge. 1) then
            call lisexp(list_load)
        end if
!
! ----- Some loads are prohibited
!
        if (l_calc) then
            call loadExcludedForAnalysis(list_load)
        end if

    end if
!
! - Modify list for undead Dirichlet loads
!
    if (l_diri_undead) then
        call load_unde_diri(list_load)
    end if

    AS_DEALLOCATE(vk8=v_list_dble)
    call jedema()
end subroutine
