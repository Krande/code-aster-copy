! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

subroutine load_neum_spec(loadName, loadNume, loadApply, ligrel_calc, i_type_neum, &
                          nb_type_neumz, nb_in_maxi, nb_in_prep, lchin, lpain, &
                          nb_in_add, load_ligrel, load_option, matr_type, iden_direct, &
                          name_inputz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
!
    character(len=8), intent(in) :: loadName
    integer, intent(in) :: loadNume
    character(len=19), intent(in) :: ligrel_calc
    character(len=4), intent(in) :: loadApply
    integer, intent(in) :: i_type_neum
    integer, intent(in) :: nb_type_neumz
    integer, intent(in) :: nb_in_maxi
    integer, intent(in) :: nb_in_prep
    character(len=*), intent(inout) :: lpain(nb_in_maxi)
    character(len=*), intent(inout) :: lchin(nb_in_maxi)
    integer, intent(out) :: nb_in_add
    character(len=19), intent(out) :: load_ligrel
    character(len=16), intent(out) :: load_option
    character(len=8), optional, intent(out) :: matr_type
    character(len=*), optional, intent(in) :: iden_direct
    character(len=*), optional, intent(in) :: name_inputz
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation
!
! Get information about load (Neumann)
!
! --------------------------------------------------------------------------------------------------
!
! In  loadName       : name of current load
! In  loadNume       : identification of load type
! In  loadApply      : type of application for load
!                        'Dead' - Dead loads (not dependent on displacements)
!                        'Pilo' - Loads for continuation (not dependent on displacements)
!                        'Suiv' - Undead loads (dependent on displacements)
! In  ligrel_calc    : LIGREL to compute
! In  i_type_neum    : index for Neumann load type
! In  nb_type_neumz  : maximum number of Neumann load type
! In  nb_in_maxi     : maximum number of input fields
! In  nb_in_prep     : number of input fields before specific ones
! IO  lpain          : list of input parameters
! IO  lchin          : list of input fields
! Out nb_in_add      : number of input fields which been added
! Out load_ligrel    : name of LIGREL for current load
! Out load_option    : name of option for current load
! Out matr_type      : matrix type for undead load
! In  iden_direct    : direct identification of type
! In  name_inputz    : direct name of input field
!
! Exception: if type VECT_ASSE -> option ='Copy_Load' and VECT_ASSE name in load_ligrel
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter(nb_type_neum=20)
    character(len=6) :: object(nb_type_neum)
    character(len=7) :: para_r(nb_type_neum), para_f(nb_type_neum)
    character(len=16) :: option_f(nb_type_neum), option_r(nb_type_neum)
    character(len=16) :: optsui_f(nb_type_neum), optsui_r(nb_type_neum)
    character(len=16) :: optmat_f(nb_type_neum), optmat_r(nb_type_neum)
    logical :: l_pilo(nb_type_neum)
    logical :: l_suiv(nb_type_neum)
    character(len=7) :: para_matr(nb_type_neum)
!
    integer :: iret, i_field_in, ig
    character(len=24) :: identify
    character(len=19) :: ligrel_load, name_input, tych, ligr_cham
    character(len=8) :: affcha, ng
    logical :: l_constant, l_fonct_0, l_fonct_t, l_sigm_int, l_find
    character(len=8), pointer :: p_vale_sigm(:) => null()
    character(len=8), pointer :: p_vect_asse(:) => null()
    integer, pointer :: desc(:) => null()
!
! - Object name construct in AFFE_CHAR_MECA
!
    data object/'.FORNO', '.F3D3D', '.F2D3D', '.F1D3D', &
        '.F2D2D', '.F1D2D', '.F1D1D', '.PESAN', &
        '.ROTAT', '.PRESS', '.FELEC', '.FCO3D', &
        '.FCO2D', '.EPSIN', '.FLUX', '.VEASS', &
        '.SIINT', '.EFOND', '.ETHM', '.ETHMH'/
!
! - Name of input parameter field (real coefficient)
!
    data para_r/'PFORNOR', 'PFR3D3D', 'PFR2D3D', 'PFR1D3D', &
        'PFR2D2D', 'PFR1D2D', 'PFR1D1D', 'PPESANR', &
        'PROTATR', 'PPRESSR', 'PFRELEC', 'PFRCO3D', &
        'PFRCO2D', 'PEPSINR', 'PFLUXR', '       ', &
        '       ', 'PEFOND', 'PECHTHM', 'HECHTHM'/
!
! - Name of option for dead load (real coefficient)
!
    data option_r/'CHAR_MECA_FORC_R', 'CHAR_MECA_FR3D3D', 'CHAR_MECA_FR2D3D', 'CHAR_MECA_FR1D3D', &
        'CHAR_MECA_FR2D2D', 'CHAR_MECA_FR1D2D', 'CHAR_MECA_FR1D1D', 'CHAR_MECA_PESA_R', &
        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_R', 'CHAR_MECA_FRELEC', 'CHAR_MECA_FRCO3D', &
        'CHAR_MECA_FRCO2D', 'CHAR_MECA_EPSI_R', 'CHAR_MECA_FLUX_R', 'Copy_Load', &
        'FORC_NODA       ', 'CHAR_MECA_EFON_R', 'CHAR_ECHA_THM_R', 'CHAR_ECHA_HR_R'/
!
! - Name of option for undead load (real coefficient)
!
    data optsui_r/'No_Load         ', 'No_Load         ', 'CHAR_MECA_FRSU23', 'No_Load         ', &
        'No_Load         ', 'No_Load         ', 'CHAR_MECA_SR1D1D', 'CHAR_MECA_PESA_R', &
        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRSU_R', 'No_Load         ', 'CHAR_MECA_SRCO3D', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'CHAR_MECA_EFSU_R', 'CHAR_ECHA_THM_R', 'CHAR_ECHA_HR_R'/
!
! - Name of option for undead load matrix (real coefficient)
!
    data optmat_r/'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'RIGI_MECA_RO    ', 'RIGI_MECA_PRSU_R', 'No_Load         ', 'RIGI_MECA_SRCO3D', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'RIGI_MECA_EFSU_R', 'No_Load         ', 'No_Load         '/
!
! - Name of input parameter field (function coefficient)
!
    data para_f/'PFORNOF', 'PFF3D3D', 'PFF2D3D', 'PFF1D3D', &
        'PFF2D2D', 'PFF1D2D', 'PFF1D1D', 'PPESANR', &
        'PROTATR', 'PPRESSF', 'PFRELEC', 'PFFCO3D', &
        'PFFCO2D', 'PEPSINF', 'PFLUXF', '       ', &
        '       ', 'PEFOND', 'PCHTHMF', 'HCHTHMF'/
!
! - Name of option for dead load (function coefficient)
!
    data option_f/'CHAR_MECA_FORC_F', 'CHAR_MECA_FF3D3D', 'CHAR_MECA_FF2D3D', 'CHAR_MECA_FF1D3D', &
        'CHAR_MECA_FF2D2D', 'CHAR_MECA_FF1D2D', 'CHAR_MECA_FF1D1D', 'CHAR_MECA_PESA_R', &
        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_F', 'CHAR_MECA_FRELEC', 'CHAR_MECA_FFCO3D', &
        'CHAR_MECA_FFCO2D', 'CHAR_MECA_EPSI_F', 'CHAR_MECA_FLUX_F', 'Copy_Load', &
        'FORC_NODA       ', 'CHAR_MECA_EFON_F', 'CHAR_ECHA_THM_F', 'CHAR_ECHA_HR_F'/
!
! - Name of option for undead load (function coefficient)
!
    data optsui_f/'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'No_Load         ', 'CHAR_MECA_SF1D1D', 'No_Load         ', &
        'No_Load         ', 'CHAR_MECA_PRSU_F', 'No_Load         ', 'CHAR_MECA_SFCO3D', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'CHAR_MECA_EFSU_F', 'CHAR_ECHA_THM_F', 'CHAR_ECHA_HR_F'/
!
! - Name of option for undead load matrix (function coefficient)
!
    data optmat_f/'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'RIGI_MECA_PRSU_F', 'No_Load         ', 'RIGI_MECA_SFCO3D', &
        'No_Load         ', 'No_Load         ', 'No_Load         ', 'No_Load         ', &
        'No_Load         ', 'RIGI_MECA_EFSU_F', 'No_Load         ', 'No_Load         '/
!
! - Flag if load can been undead load type
!
    data l_suiv/.false., .false., .true., .false., &
        .false., .false., .true., .true., &
        .true., .true., .false., .true., &
        .false., .false., .false., .false., &
        .false., .true., .true., .true./
!
! - Flag if load can been used for continuation methods
!
    data l_pilo/.true., .true., .true., .true., &
        .true., .true., .true., .true., &
        .false., .true., .false., .true., &
        .true., .false., .true., .true., &
        .false., .true., .false., .false./
!
! - Type of matrix for undead load
!
    data para_matr/'       ', '       ', '       ', '       ', &
        '       ', '       ', '       ', '       ', &
        'PMATUUR', 'PMATUNS', '       ', 'PMATUNS', &
        '       ', '       ', '       ', '       ', &
        '       ', 'PMATUNS', '       ', '       '/
!
! --------------------------------------------------------------------------------------------------
!
    ligrel_load = loadName(1:8)//'.CHME.LIGRE'
    load_ligrel = ' '
    load_option = 'No_Load'
    l_constant = .false.
    l_fonct_0 = .false.
    l_fonct_t = .false.
    l_sigm_int = .false.
    i_field_in = nb_in_prep
    ASSERT(i_type_neum .le. nb_type_neum)
    ASSERT(nb_type_neumz .eq. nb_type_neum)
!
! - Identify current load
!
    iret = 0
    if (present(iden_direct)) then
        if (iden_direct .eq. object(i_type_neum)) then
            iret = 1
        end if
        name_input = name_inputz
    else
        if (object(i_type_neum) .eq. '.VEASS') then
            identify = loadName(1:8)//'.CHME'//object(i_type_neum)
        else
            identify = loadName(1:8)//'.CHME'//object(i_type_neum)//'.DESC'
        end if
        name_input = loadName(1:8)//'.CHME'//object(i_type_neum)
        call jeexin(identify, iret)
    end if
!
    if (iret .ne. 0) then
!
! ----- Value type
!
        if (loadNume .eq. 1) then
            l_constant = .true.
        elseif (loadNume .eq. 5) then
            l_constant = .true.
        elseif (loadNume .eq. 8) then
            l_fonct_0 = .true.
        else if (loadNume .eq. 2) then
            l_fonct_0 = .true.
        else if (loadNume .eq. 3) then
            l_fonct_t = .true.
        else if (loadNume .eq. 55) then
            l_sigm_int = .true.
        end if
!
! ----- Special for undeads loads
!
        if (loadNume .eq. 4 .or. loadNume .eq. 9 .or. loadNume .eq. 11) then
            l_constant = .true.
            call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=affcha)
            if (affcha(5:7) .eq. '_FO') then
                l_constant = .false.
                l_fonct_0 = .true.
            end if
        end if
!
! ----- Name of option
!
        if (l_constant) then
            if (loadApply .eq. 'Suiv') then
                if (present(matr_type)) then
                    load_option = optmat_r(i_type_neum)
                else
                    load_option = optsui_r(i_type_neum)
                end if
            else
                load_option = option_r(i_type_neum)
            end if
        else if (l_fonct_0 .or. l_fonct_t) then
            if (loadApply .eq. 'Suiv') then
                if (present(matr_type)) then
                    load_option = optmat_f(i_type_neum)
                else
                    load_option = optsui_f(i_type_neum)
                end if
            else
                load_option = option_f(i_type_neum)
            end if
        else if (l_sigm_int) then
            load_option = option_r(i_type_neum)
            ASSERT(load_option .eq. 'FORC_NODA')
        else
            ASSERT(.false.)
        end if
!
! ----- Name of input fields
!
        if (l_constant) then
            i_field_in = i_field_in+1
            lpain(i_field_in) = para_r(i_type_neum)
            lchin(i_field_in) = name_input(1:19)
!
            if (lpain(i_field_in) .eq. 'PEPSINR') then
                call jeveuo(ligrel_load(1:13)//'.EPSIN.DESC', 'L', vi=desc)
                ig = desc(1)
                call jenuno(jexnum('&CATA.GD.NOMGD', ig), ng)
!               recuperation du nom du champ stocké dans la carte "bidon"
                if (ng .eq. 'NEUT_K8') then
                    call jeveuo(ligrel_load(1:13)//'.EPSIN.VALE', 'L', vk8=p_vale_sigm)
                    lchin(i_field_in) = p_vale_sigm(1)
                end if
            end if

            if (load_option .eq. 'CHAR_MECA_EFON_R') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFR'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            else if (load_option .eq. 'CHAR_MECA_EFSU_R') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFR'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            else if (load_option .eq. 'RIGI_MECA_EFSU_R') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFR'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            else if (load_option .eq. 'CHAR_ECHA_THM_R') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PDEPLMR'
                lchin(i_field_in) = loadName//'.CHME.DEPL_R'
            else if (load_option .eq. 'CHAR_ECHA_HR_R') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PDEPLMR'
                lchin(i_field_in) = loadName//'.CHME.DEPL_R'
            else if (load_option .eq. 'CHAR_ECHA_THM_F') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PDEPLMR'
                lchin(i_field_in) = loadName//'.CHME.DEPL_R'
            else if (load_option .eq. 'CHAR_ECHA_HR_F') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PDEPLMR'
                lchin(i_field_in) = loadName//'.CHME.DEPL_R'
            end if
        else if (l_fonct_0 .or. l_fonct_t) then
            i_field_in = i_field_in+1
            lpain(i_field_in) = para_f(i_type_neum)
            lchin(i_field_in) = name_input(1:19)
            if (load_option .eq. 'CHAR_MECA_EFON_F') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFF'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            end if
            if (load_option .eq. 'CHAR_MECA_EFSU_F') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFF'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            end if
            if (load_option .eq. 'RIGI_MECA_EFSU_F') then
                i_field_in = i_field_in+1
                lpain(i_field_in) = 'PPREFFF'
                lchin(i_field_in) = loadName//'.CHME.PREFF'
            end if
        else if (l_sigm_int) then
            call jeveuo(ligrel_load(1:13)//'.SIINT.VALE', 'L', vk8=p_vale_sigm)
            i_field_in = i_field_in+1
            lpain(i_field_in) = 'PSIEFR'
            lchin(i_field_in) = p_vale_sigm(1)
        else
            ASSERT(.false.)
        end if
!
! ----- Parameter name of output field for matrix (undead loads)
!
        if (loadApply .eq. 'Suiv') then
            if (present(matr_type)) then
                matr_type = para_matr(i_type_neum)
            end if
        end if
!
! ----- Select LIGREL
!
        if (object(i_type_neum) .eq. '.FORNO') then
            load_ligrel = ligrel_load
        else if (present(name_inputz)) then
            l_find = ASTER_FALSE
            call dismoi('TYPE_CHAMP', name_input, 'CHAMP', repk=tych)
            if (tych(1:2) == "EL") then
                call dismoi('NOM_LIGREL', name_input, 'CHAMP', repk=ligr_cham)
                call exisd('LIGREL', ligr_cham, iret)
                if (iret == 1) then
                    load_ligrel = ligr_cham
                    l_find = ASTER_TRUE
                end if
            end if
            if (.not. l_find) then
                load_ligrel = ligrel_calc
            end if
        else
            load_ligrel = ligrel_calc
        end if
        if (load_option .eq. 'Copy_Load') then
            ASSERT((loadNume .ge. 1 .and. loadNume .le. 3) .or. loadNume .eq. 5)
            call jeveuo(lchin(i_field_in), 'L', vk8=p_vect_asse)
            load_ligrel = p_vect_asse(1)
        end if
!
! ----- Checking for undead loads
!
        if (loadApply .eq. 'Suiv') then
            if (.not. l_suiv(i_type_neum)) then
                call utmess('F', 'CHARGES_23', sk=loadName)
            end if
            if ((load_option .eq. 'No_Load') .and. (.not. present(matr_type))) then
                call utmess('F', 'CHARGES_23', sk=loadName)
            end if
        end if
!
! ----- Checking for continuation type loads
!
        if (loadApply .eq. 'Pilo') then
            if (.not. l_pilo(i_type_neum)) then
                call utmess('F', 'CHARGES_26', sk=loadName)
            end if
            if (l_fonct_t) then
                call utmess('F', 'CHARGES_28', sk=loadName)
            end if
        end if
!
! ----- Number of input fields which been added
!
        nb_in_add = i_field_in-nb_in_prep
!
        ASSERT(i_field_in .le. nb_in_maxi)
    end if
!
end subroutine
