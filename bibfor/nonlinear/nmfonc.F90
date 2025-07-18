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
subroutine nmfonc(ds_conv, ds_algopara, solver, model, ds_contact, &
                  list_load, sdnume, sddyna, ds_errorindic, mater, &
                  ds_inout, ds_constitutive, ds_energy, ds_algorom, ds_posttimestep, &
                  list_func_acti)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/GetResi.h"
#include "asterfort/infniv.h"
#include "asterfort/ischar.h"
#include "asterfort/isdiri.h"
#include "asterfort/isfonc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/matdis.h"
#include "asterfort/ndynlo.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_Conv), intent(in) :: ds_conv
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19), intent(in) :: solver
    character(len=24), intent(in) :: model
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: list_load
    character(len=19), intent(in) :: sdnume
    character(len=19), intent(in) :: sddyna
    type(NL_DS_ErrorIndic), intent(in) :: ds_errorindic
    character(len=24), intent(in) :: mater
    type(NL_DS_InOut), intent(in) :: ds_inout
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Energy), intent(in) :: ds_energy
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_PostTimeStep), intent(in) :: ds_posttimestep
    integer(kind=8), intent(inout) :: list_func_acti(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Prepare active functionnalities information
!
! --------------------------------------------------------------------------------------------------
!
! NB: to ask list_func_acti, use ISFONC.F90 subroutine !
!
! In  ds_conv          : datastructure for convergence management
! In  ds_algopara      : datastructure for algorithm parameters
! In  solver           : datastructure for solver parameters
! In  model            : name of the model
! In  ds_contact       : datastructure for contact management
! In  list_load        : name of datastructure for list of loads
! In  sdnume           : datastructure for dof positions
! In  sddyna           : name of datastructure for dynamic parameters
! In  ds_errorindic    : datastructure for error indicator
! In  mater            : name of material
! In  ds_inout         : datastructure for input/output management
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_energy        : datastructure for energy management
! In  ds_algorom       : datastructure for ROM parameters
! In  ds_posttimestep  : datastructure for post-treatment at each time step
! IO  list_func_acti   : list of active functionnalities
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nocc, iret, nb_subs_stat, nb_sst
    integer(kind=8) :: i_cont_form
    aster_logical :: l_deborst, l_frot, l_dis_choc, l_all_verif, l_refe, l_comp, lAnnealing
    aster_logical :: l_loop_geom, l_loop_frot, l_loop_cont, l_pena, l_rela, l_maxi
    integer(kind=8) :: ixfem
    aster_logical :: l_load_undead, l_load_elim, l_load_didi
    character(len=8) :: k8bid, repk
    character(len=16) :: command, k16bid
    character(len=3)  :: matr_distr
    character(len=24) :: solv_type, solv_precond
    aster_logical :: l_stat, l_dyna
    aster_logical :: l_newt_cont, l_newt_frot, l_newt_geom
    aster_logical :: l_dyna_expl, l_cont, l_unil
    character(len=24), pointer :: slvk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    l_cont = ds_contact%l_meca_cont
    l_unil = ds_contact%l_meca_unil
    l_deborst = ds_constitutive%l_deborst
    l_dis_choc = ds_constitutive%l_dis_choc
    lAnnealing = ds_constitutive%lAnnealing
!
! - Print
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_11')
    end if
!
! - Command
!
    call getres(k8bid, k16bid, command)
    l_stat = command(1:4) .eq. 'STAT'
    l_dyna = command(1:4) .eq. 'DYNA'
    l_dyna_expl = ndynlo(sddyna, 'EXPLICITE')
!
! - Large rotations
!
    call jeexin(sdnume(1:19)//'.NDRO', iret)
    if (iret .gt. 0) list_func_acti(15) = 1
!
! - Damaged nodes
!
    call jeexin(sdnume(1:19)//'.ENDO', iret)
    if (iret .gt. 0) list_func_acti(40) = 1
!
! - Line search
!
    if (ds_algopara%l_line_search) then
        list_func_acti(1) = 1
    end if
!
! - Continuation methods (PILOTAGE)
!
    if (l_stat) then
        call getfac('PILOTAGE', nocc)
        if (nocc .ne. 0) list_func_acti(2) = 1
    end if
!
! - Unilateral condition
!
    if (l_unil) list_func_acti(12) = 1
!
! - Energy computation
!
    if (ds_energy%l_comp) then
        list_func_acti(50) = 1
    end if
!
! - Dynamic
!
    if (l_dyna) then
        list_func_acti(69) = 1
    end if
!
! - Modal projection for dynamic
!
    if (l_dyna) then
        if (ndynlo(sddyna, 'PROJ_MODAL')) then
            list_func_acti(51) = 1
        end if
    end if
!
! - Distributed matrix (parallel computaing)
!
    call matdis(matr_distr)
    if (matr_distr .eq. 'OUI') list_func_acti(52) = 1
!
! - Deborst algorithm
!
    if (l_deborst) then
        list_func_acti(7) = 1
    end if
!
! - Reference criterion RESI_REFE_RELA
!
    call GetResi(ds_conv, type='RESI_REFE_RELA', l_resi_test_=l_refe)
    if (l_refe) then
        list_func_acti(8) = 1
    end if
!
! - By components criterion RESI_COMP_RELA
!
    call GetResi(ds_conv, type='RESI_COMP_RELA', l_resi_test_=l_comp)
    if (l_comp) then
        list_func_acti(35) = 1
    end if
!
! - Global relative criterion RESI_GLOB_RELA
!
    call GetResi(ds_conv, type='RESI_GLOB_RELA', l_resi_test_=l_rela)
    if (l_rela) then
        list_func_acti(70) = 1
    end if
!
! - Global maximal criterion RESI_GLOB_MAXI
!
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', l_resi_test_=l_maxi)
    if (l_maxi) then
        list_func_acti(71) = 1
    end if
!
! - X-FEM
!
    call exixfe(model, ixfem)
    if (ixfem .ne. 0) list_func_acti(6) = 1
!
! - Contact_friction
!
    if (l_cont) then
        i_cont_form = cfdisi(ds_contact%sdcont_defi, 'FORMULATION')
        list_func_acti(64) = 1
        if (i_cont_form .eq. 2) then
            list_func_acti(5) = 1
            list_func_acti(17) = cfdisi(ds_contact%sdcont_defi, 'ALL_INTERPENETRE')
            list_func_acti(26) = 1
            l_frot = cfdisl(ds_contact%sdcont_defi, 'FROTTEMENT')
            if (l_frot) then
                list_func_acti(10) = 1
                list_func_acti(27) = 1
            end if
        else if (i_cont_form .eq. 1) then
            list_func_acti(4) = 1
            l_frot = cfdisl(ds_contact%sdcont_defi, 'FROTTEMENT')
            if (l_frot) then
                list_func_acti(3) = 1
            end if
        else if (i_cont_form .eq. 5) then
            list_func_acti(63) = 1
            list_func_acti(26) = 1
            l_frot = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
! - Contact: no computation
!
    if (l_cont) then
        l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
        if (l_all_verif) list_func_acti(38) = 1
    end if
!
! - Contact: fixed loops
!
    if (l_cont) then
        l_loop_geom = ASTER_FALSE
        l_loop_frot = ASTER_FALSE
        l_loop_cont = ASTER_FALSE
        l_loop_geom = cfdisl(ds_contact%sdcont_defi, 'GEOM_BOUCLE')
        if (l_frot) l_loop_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_BOUCLE')
        l_loop_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_BOUCLE')
        if (l_all_verif) then
            l_loop_cont = ASTER_FALSE
            l_loop_geom = ASTER_FALSE
            l_loop_frot = ASTER_FALSE
        end if
        if (i_cont_form .eq. 1) then
            l_loop_cont = ASTER_FALSE
            l_loop_frot = ASTER_FALSE
        end if
        if (l_loop_geom) list_func_acti(31) = 1
        if (l_loop_frot) list_func_acti(32) = 1
        if (l_loop_cont) list_func_acti(33) = 1
        if (l_loop_geom .or. l_loop_frot .or. l_loop_cont) then
            list_func_acti(34) = 1
        end if
    end if
!
! - Generalized Newton
!
    if (l_cont) then
        if (i_cont_form .eq. 2 .or. i_cont_form .eq. 5) then
            l_newt_geom = cfdisl(ds_contact%sdcont_defi, 'GEOM_NEWTON')
            l_newt_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_NEWTON')
            l_newt_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_NEWTON')
            if (l_newt_frot) list_func_acti(47) = 1
            if (l_newt_cont) list_func_acti(53) = 1
            if (l_newt_geom) list_func_acti(55) = 1
        end if
    end if
! At least one contact zone has penalisation method
    if (l_cont) then
        if (i_cont_form .eq. 2 .or. i_cont_form .eq. 5) then
            l_pena = cfdisl(ds_contact%sdcont_defi, 'EXIS_PENA')
            if (l_pena) list_func_acti(66) = 1
        end if
    end if
!
! - At least, one Neumann undead load ?
!
    l_load_undead = ischar(list_load, 'NEUM', 'SUIV')
    if (l_load_undead) then
        list_func_acti(13) = 1
    end if
!
! - At least, one Dirichlet undead load ?
!
    l_load_undead = ischar(list_load, 'DIRI', 'SUIV')
    if (l_load_undead) then
        list_func_acti(60) = 1
    end if
!
! - At least, one "DIDI" load ?
!
    l_load_didi = ischar(list_load, 'DIRI', 'DIDI')
    if (l_load_didi) then
        list_func_acti(22) = 1
    end if
!
! - At least, one AFFE_CHAR_CINE load ?
!
    l_load_elim = isdiri(list_load, 'ELIM')
    if (l_load_elim) then
        list_func_acti(36) = 1
    end if
!
! - Static substructuring
!
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nb_subs_stat)
    if (nb_subs_stat .gt. 0) list_func_acti(14) = 1
!
! - Substructuring loads
!
    call getfac('SOUS_STRUC', nb_sst)
    if (nb_sst .gt. 0) then
        list_func_acti(24) = 1
    end if
!
! - Buckling
!
    if (ds_posttimestep%l_crit_stab) then
        list_func_acti(18) = 1
    end if
!
! - Stability
!
    if (ds_posttimestep%stab_para%nb_dof_stab .gt. 0) then
        list_func_acti(49) = 1
    end if
!
! - Vibration modes
!
    if (ds_posttimestep%l_mode_vibr) then
        list_func_acti(19) = 1
    end if
!
! - THM time error
!
    if (l_stat) then
        if (ds_errorindic%l_erre_thm) then
            list_func_acti(21) = 1
        end if
    end if
!
! - IMPLEX algorithm
!
    if (l_stat) then
        if (ds_algopara%method .eq. 'IMPLEX') then
            list_func_acti(28) = 1
        end if
    end if
!
! - NEWTON_KRYLOV algorithm
!
    if (ds_algopara%method .eq. 'NEWTON_KRYLOV') then
        list_func_acti(48) = 1
    end if
!
! - NEWTON_KRYLOV algorithm
!
    if (ds_algopara%method .eq. 'MODELE_REDUIT') then
        list_func_acti(61) = 1
        if (ds_algorom%l_hrom) then
            list_func_acti(62) = 1
            if (ds_algorom%l_hrom_corref) then
                list_func_acti(67) = 1
                list_func_acti(34) = 1
            end if
        end if
    end if
!
! - DIS_CHOC elements ?
!
    if (l_dis_choc) then
        list_func_acti(29) = 1
    end if
!
! - Command variables
!
    call dismoi('EXI_VARC', mater, 'CHAM_MATER', repk=repk)
    if (repk .eq. 'OUI') list_func_acti(30) = 1
!
! - THM ?
!
    call dismoi('EXI_THM', model, 'MODELE', repk=repk)
    if (repk .eq. 'OUI') then
        list_func_acti(37) = 1
    end if

!
! - Elemesnt with STRX field (multifibers for instantce)
!
    call dismoi('EXI_STRX', model, 'MODELE', repk=repk)
    if (repk .eq. 'OUI') list_func_acti(56) = 1
!
! - REUSE ?
!
    if (ds_inout%l_reuse) then
        list_func_acti(39) = 1
    end if
!
! - Does ETAT_INIT (initial state) exist ?
!
    if (ds_inout%l_state_init) then
        list_func_acti(59) = 1
    end if
!
! - Solvers
!
    call jeveuo(solver//'.SLVK', 'L', vk24=slvk)
    solv_type = slvk(1)
    if (solv_type .eq. 'LDLT') list_func_acti(41) = 1
    if (solv_type .eq. 'MULT_FRONT') list_func_acti(42) = 1
    if (solv_type .eq. 'GCPC') list_func_acti(43) = 1
    if (solv_type .eq. 'MUMPS') list_func_acti(44) = 1
    if (solv_type .eq. 'PETSC') list_func_acti(45) = 1
    if (solv_type .eq. 'PETSC' .or. solv_type .eq. 'GCPC') then
        solv_precond = slvk(2)
        if (solv_precond .eq. 'LDLT_SP') list_func_acti(46) = 1
    end if
!
! - Explicit dynamics
!
    if (l_dyna_expl) list_func_acti(54) = 1
!
! - Do elastic properties are functions ?
!
    call dismoi('ELAS_FO', mater, 'CHAM_MATER', repk=repk)
    if (repk .eq. 'OUI') list_func_acti(57) = 1

! - Annealing ?
    if (lAnnealing) then
        list_func_acti(58) = 1
    end if
!
! - HHO ?
!
    call dismoi('EXI_HHO', model, 'MODELE', repk=repk)
    if (repk .eq. 'OUI') then
        list_func_acti(68) = 1
    end if
!
! - Print
!
    if (niv .ge. 2) then
!
! ----- Solving methods
!
        if (isfonc(list_func_acti, 'IMPLEX')) then
            call utmess('I', 'MECANONLINE14_1')
        end if
        if (isfonc(list_func_acti, 'EXPLICITE')) then
            call utmess('I', 'MECANONLINE14_2')
        end if
        if (isfonc(list_func_acti, 'DYNAMIQUE')) then
            call utmess('I', 'MECANONLINE14_65')
        end if
        if (isfonc(list_func_acti, 'NEWTON_KRYLOV')) then
            call utmess('I', 'MECANONLINE14_3')
        end if
        if (isfonc(list_func_acti, 'ROM')) then
            call utmess('I', 'MECANONLINE14_4')
        end if
        if (isfonc(list_func_acti, 'HROM')) then
            call utmess('I', 'MECANONLINE14_5')
        end if
        if (isfonc(list_func_acti, 'HROM_CORR_EF')) then
            call utmess('I', 'MECANONLINE14_6')
        end if
        if (isfonc(list_func_acti, 'RECH_LINE')) then
            call utmess('I', 'MECANONLINE14_7')
        end if
        if (isfonc(list_func_acti, 'PILOTAGE')) then
            call utmess('I', 'MECANONLINE14_8')
        end if
        if (isfonc(list_func_acti, 'DEBORST')) then
            call utmess('I', 'MECANONLINE14_9')
        end if
        if (isfonc(list_func_acti, 'SOUS_STRUC')) then
            call utmess('I', 'MECANONLINE14_10')
        end if
        if (isfonc(list_func_acti, 'PROJ_MODAL')) then
            call utmess('I', 'MECANONLINE14_11')
        end if
!
! ----- Contact
!
        if (isfonc(list_func_acti, 'CONTACT')) then
            call utmess('I', 'MECANONLINE14_12')
        end if
        if (isfonc(list_func_acti, 'CONT_DISCRET')) then
            call utmess('I', 'MECANONLINE14_13')
        end if
        if (isfonc(list_func_acti, 'CONT_CONTINU')) then
            call utmess('I', 'MECANONLINE14_14')
        end if
        if (isfonc(list_func_acti, 'CONT_LAC')) then
            call utmess('I', 'MECANONLINE14_17')
        end if
        if (isfonc(list_func_acti, 'BOUCLE_EXT_GEOM')) then
            call utmess('I', 'MECANONLINE14_18')
        end if
        if (isfonc(list_func_acti, 'BOUCLE_EXT_CONT')) then
            call utmess('I', 'MECANONLINE14_19')
        end if
        if (isfonc(list_func_acti, 'BOUCLE_EXT_FROT')) then
            call utmess('I', 'MECANONLINE14_20')
        end if
        if (isfonc(list_func_acti, 'BOUCLE_EXTERNE')) then
            call utmess('I', 'MECANONLINE14_21')
        end if
        if (isfonc(list_func_acti, 'GEOM_NEWTON')) then
            call utmess('I', 'MECANONLINE14_22')
        end if
        if (isfonc(list_func_acti, 'FROT_NEWTON')) then
            call utmess('I', 'MECANONLINE14_23')
        end if
        if (isfonc(list_func_acti, 'CONT_NEWTON')) then
            call utmess('I', 'MECANONLINE14_24')
        end if
        if (isfonc(list_func_acti, 'CONT_ALL_VERIF')) then
            call utmess('I', 'MECANONLINE14_25')
        end if
        if (isfonc(list_func_acti, 'CONTACT_INIT')) then
            call utmess('I', 'MECANONLINE14_26')
        end if
        if (isfonc(list_func_acti, 'LIAISON_UNILATER')) then
            call utmess('I', 'MECANONLINE14_27')
        end if
        if (isfonc(list_func_acti, 'FROT_DISCRET')) then
            call utmess('I', 'MECANONLINE14_28')
        end if
        if (isfonc(list_func_acti, 'FROT_CONTINU')) then
            call utmess('I', 'MECANONLINE14_29')
        end if
        if (isfonc(list_func_acti, 'EXIS_PENA')) then
            call utmess('I', 'MECANONLINE14_64')
        end if
!
! ----- Finite elements
!
        if (isfonc(list_func_acti, 'ELT_CONTACT')) then
            call utmess('I', 'MECANONLINE14_31')
        end if
        if (isfonc(list_func_acti, 'DIS_CHOC')) then
            call utmess('I', 'MECANONLINE14_33')
        end if
        if (isfonc(list_func_acti, 'GD_ROTA')) then
            call utmess('I', 'MECANONLINE14_34')
        end if
        if (isfonc(list_func_acti, 'XFEM')) then
            call utmess('I', 'MECANONLINE14_35')
        end if
        if (isfonc(list_func_acti, 'EXI_STRX')) then
            call utmess('I', 'MECANONLINE14_36')
        end if
!
! ----- CONVERGENCE
!
        if (isfonc(list_func_acti, 'RESI_REFE')) then
            call utmess('I', 'MECANONLINE14_37')
        end if
        if (isfonc(list_func_acti, 'RESI_COMP')) then
            call utmess('I', 'MECANONLINE14_38')
        end if
!
! ----- Loads
!
        if (isfonc(list_func_acti, 'NEUM_UNDEAD')) then
            call utmess('I', 'MECANONLINE14_39')
        end if
        if (isfonc(list_func_acti, 'DIRI_UNDEAD')) then
            call utmess('I', 'MECANONLINE14_40')
        end if
        if (isfonc(list_func_acti, 'DIDI')) then
            call utmess('I', 'MECANONLINE14_41')
        end if
        if (isfonc(list_func_acti, 'DIRI_CINE')) then
            call utmess('I', 'MECANONLINE14_42')
        end if
!
! ----- MODELISATION
!
        if (isfonc(list_func_acti, 'MACR_ELEM_STAT')) then
            call utmess('I', 'MECANONLINE14_44')
        end if
        if (isfonc(list_func_acti, 'THM')) then
            call utmess('I', 'MECANONLINE14_45')
        end if
        if (isfonc(list_func_acti, 'HHO')) then
            call utmess('I', 'MECANONLINE14_32')
        end if
        if (isfonc(list_func_acti, 'ENDO_NO')) then
            call utmess('I', 'MECANONLINE14_46')
        end if
!
! ----- Post-treatments
!
        if (isfonc(list_func_acti, 'CRIT_STAB')) then
            call utmess('I', 'MECANONLINE14_47')
        end if
        if (isfonc(list_func_acti, 'DDL_STAB')) then
            call utmess('I', 'MECANONLINE14_48')
        end if
        if (isfonc(list_func_acti, 'MODE_VIBR')) then
            call utmess('I', 'MECANONLINE14_49')
        end if
        if (isfonc(list_func_acti, 'ENERGIE')) then
            call utmess('I', 'MECANONLINE14_50')
        end if
        if (isfonc(list_func_acti, 'ERRE_TEMPS_THM')) then
            call utmess('I', 'MECANONLINE14_51')
        end if
        if (isfonc(list_func_acti, 'ANNEALING')) then
            call utmess('I', 'MECANONLINE14_52')
        end if
        if (isfonc(list_func_acti, 'EXI_VARC')) then
            call utmess('I', 'MECANONLINE14_53')
        end if
        if (isfonc(list_func_acti, 'ELAS_FO')) then
            call utmess('I', 'MECANONLINE14_54')
        end if
!
        if (isfonc(list_func_acti, 'REUSE')) then
            call utmess('I', 'MECANONLINE14_55')
        end if
        if (isfonc(list_func_acti, 'ETAT_INIT')) then
            call utmess('I', 'MECANONLINE14_56')
        end if
!
! ----- Solver options
!
        if (isfonc(list_func_acti, 'LDLT')) then
            call utmess('I', 'MECANONLINE14_57')
        end if
        if (isfonc(list_func_acti, 'MULT_FRONT')) then
            call utmess('I', 'MECANONLINE14_58')
        end if
        if (isfonc(list_func_acti, 'GCPC')) then
            call utmess('I', 'MECANONLINE14_59')
        end if
        if (isfonc(list_func_acti, 'MUMPS')) then
            call utmess('I', 'MECANONLINE14_60')
        end if
        if (isfonc(list_func_acti, 'PETSC')) then
            call utmess('I', 'MECANONLINE14_61')
        end if
        if (isfonc(list_func_acti, 'LDLT_SP')) then
            call utmess('I', 'MECANONLINE14_62')
        end if
        if (isfonc(list_func_acti, 'MATR_DISTRIBUEE')) then
            call utmess('I', 'MECANONLINE14_63')
        end if
    end if
!
end subroutine
