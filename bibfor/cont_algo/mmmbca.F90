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
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
subroutine mmmbca(mesh, iter_newt, nume_inst, &
                  sddisc, disp_curr, disp_cumu_inst, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfdisr.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnumm.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcomce.h"
#include "asterfort/mm_cycl_algo.h"
#include "asterfort/mm_cycl_prop.h"
#include "asterfort/mm_cycl_stat.h"
#include "asterfort/mmbouc.h"
#include "asterfort/mmeval_prep.h"
#include "asterfort/mmeven.h"
#include "asterfort/mmextm.h"
#include "asterfort/mmfield_prep.h"
#include "asterfort/mmglis.h"
#include "asterfort/mmimp4.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
#include "asterfort/mminfr.h"
#include "asterfort/mmstac.h"
#include "asterfort/mmstaf.h"
#include "asterfort/mreacg.h"
#include "Contact_type.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: iter_newt
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: disp_curr
    character(len=19), intent(in) :: disp_cumu_inst
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Management of contact loop
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  iter_newt        : index of current Newton iteration
! In  nume_inst        : index of current time step
! In  sddisc           : datastructure for time discretization
! In  disp_curr        : current displacements
! In  disp_cumu_inst   : displacement increment from beginning of current time
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ztabf
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: jdecme, elem_slav_indx, elem_slav_nume, elem_mast_nume
    integer(kind=8) :: indi_cont_curr, indi_cont_prev, indi_frot_prev, indi_frot_curr
    integer(kind=8) :: i_zone, i_elem_slav, i_cont_poin, i_poin_elem
    integer(kind=8) :: model_ndim, nb_cont_zone, loop_cont_vali
    integer(kind=8) :: elem_slav_nbno, nb_poin_elem, nb_elem_slav
    integer(kind=8) :: indi_cont_eval, indi_frot_eval
    real(kind=8) :: ksipr1, ksipr2, ksipc1, ksipc2
    real(kind=8) :: ksipr1_old, ksipr2_old, ksipc1_old, ksipc2_old, resi_geom
    real(kind=8) :: norm(3), tau1(3), tau2(3)
    real(kind=8) :: lagr_cont_node(9), lagr_fro1_node(9), lagr_fro2_node(9)
    real(kind=8) :: elem_slav_coor(27)
    real(kind=8) :: lagr_cont_poin, time_curr
    real(kind=8) :: gap, gap_user
    real(kind=8) :: pres_frot(3), gap_user_frot(3)
    real(kind=8) :: coef_cont, coef_frot, loop_cont_vale
    character(len=8) :: elem_slav_type
    character(len=19) :: cnscon, cnsfr1, cnsfr2
    character(len=19) :: oldgeo, newgeo
    character(len=19) :: chdepd
    aster_logical :: l_glis
    aster_logical :: l_glis_init, l_veri, l_exis_glis, loop_cont_conv, l_loop_cont
    aster_logical :: l_frot_zone, l_pena_frot, l_frot
    aster_logical :: l_pena_cont
    integer(kind=8) :: loop_geom_count, loop_fric_count, loop_cont_count
    integer(kind=8) :: type_adap
    character(len=24) :: sdcont_cychis, sdcont_cyccoe, sdcont_cyceta
    real(kind=8), pointer :: v_sdcont_cychis(:) => null()
    real(kind=8), pointer :: v_sdcont_cyccoe(:) => null()
    integer(kind=8), pointer :: v_sdcont_cyceta(:) => null()
    character(len=24) :: sdcont_tabfin, sdcont_jsupco, sdcont_apjeu
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    real(kind=8), pointer :: v_sdcont_jsupco(:) => null()
    real(kind=8), pointer :: v_sdcont_apjeu(:) => null()
    real(kind=8)  :: vale_pene, glis_maxi, resi_cont
    real(kind=8)  :: sum_cont_press, resi_press_curr
    real(kind=8)  :: coor_escl_curr(3), coor_proj_curr(3)
    real(kind=8)  :: wpg_old
    aster_logical :: l_coef_adap
    character(len=8) :: iptxt
    integer(kind=8) :: hist_index, coun_bcle_geom, nb_cont_poin
    aster_logical :: l_granglis
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... ACTIVATION/DESACTIVATION'
    end if
!
! - Initializations
!
    l_frot_zone = ASTER_FALSE
    l_pena_frot = ASTER_FALSE
    l_frot = ASTER_FALSE
    l_pena_cont = ASTER_FALSE
    loop_cont_conv = ASTER_TRUE
    loop_cont_vali = 0
    ds_contact%critere_geom = 0.0
    resi_geom = 0.0
    vale_pene = 0.0
    glis_maxi = 0.d0
    resi_cont = -1.0
    sum_cont_press = -1.0
    resi_press_curr = 1.D6
    coor_escl_curr(:) = 0.0
    coor_proj_curr(:) = 0.0
!
! - Parameters
!
    l_exis_glis = cfdisl(ds_contact%sdcont_defi, 'EXIS_GLISSIERE')
    l_loop_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_BOUCLE')
    type_adap = cfdisi(ds_contact%sdcont_defi, 'TYPE_ADAPT')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    nb_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROTTEMENT')
    resi_cont = cfdisr(ds_contact%sdcont_defi, 'RESI_CONT')
!
! - Acces to contact objects
!
    ztabf = cfmmvd('ZTABF')
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    sdcont_jsupco = ds_contact%sdcont_solv(1:14)//'.JSUPCO'
    sdcont_apjeu = ds_contact%sdcont_solv(1:14)//'.APJEU'
    call jeveuo(sdcont_tabfin, 'E', vr=v_sdcont_tabfin)
    call jeveuo(sdcont_jsupco, 'E', vr=v_sdcont_jsupco)
    call jeveuo(sdcont_apjeu, 'E', vr=v_sdcont_apjeu)
!
! - Acces to cycling objects
!
    sdcont_cyceta = ds_contact%sdcont_solv(1:14)//'.CYCETA'
    sdcont_cychis = ds_contact%sdcont_solv(1:14)//'.CYCHIS'
    sdcont_cyccoe = ds_contact%sdcont_solv(1:14)//'.CYCCOE'
    call jeveuo(sdcont_cyceta, 'L', vi=v_sdcont_cyceta)
    call jeveuo(sdcont_cychis, 'E', vr=v_sdcont_cychis)
    call jeveuo(sdcont_cyccoe, 'E', vr=v_sdcont_cyccoe)
!
! - Get current time
!
    time_curr = diinst(sddisc, nume_inst)
    ds_contact%time_curr = time_curr
!
! - Geometric update
!
    oldgeo = mesh//'.COORDO'
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
    call mreacg(mesh, ds_contact, disp_curr)
!
! - Prepare displacement field to get contact Lagrangien multiplier
!
    cnscon = '&&MMMBCA.CNSCON'
    call mmfield_prep(disp_curr, cnscon, &
                      l_sort_=ASTER_TRUE, nb_cmp_=1, list_cmp_=['LAGS_C  '])
!
! - Prepare displacement field to get friction Lagrangien multiplier
!
    chdepd = '&&MMMBCA.CHDEPD'
    cnsfr1 = '&&MMMBCA.CNSFR1'
    cnsfr2 = '&&MMMBCA.CNSFR2'
    if (l_frot) then
        call mmfield_prep(disp_cumu_inst, cnsfr1, &
                          l_sort_=ASTER_TRUE, nb_cmp_=1, list_cmp_=['LAGS_F1 '])
        if (model_ndim .eq. 3) then
            call mmfield_prep(disp_cumu_inst, cnsfr2, &
                              l_sort_=ASTER_TRUE, nb_cmp_=1, list_cmp_=['LAGS_F2 '])
        end if
        call mmfield_prep(oldgeo, chdepd, &
                          l_update_=ASTER_TRUE, field_update_=disp_cumu_inst)
    end if
!
! - Loop on contact zones
!
    sum_cont_press = 0.d0
    ds_contact%resi_pressure = 0.d0
    i_cont_poin = 1
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters of zone
!
        l_glis = mminfl(ds_contact%sdcont_defi, 'GLISSIERE_ZONE', i_zone)
        l_veri = mminfl(ds_contact%sdcont_defi, 'VERIF', i_zone)
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', i_zone)
        l_frot_zone = mminfl(ds_contact%sdcont_defi, 'FROTTEMENT_ZONE', i_zone)
        l_pena_frot = mminfl(ds_contact%sdcont_defi, 'ALGO_FROT_PENA', i_zone)
        l_pena_cont = mminfl(ds_contact%sdcont_defi, 'ALGO_CONT_PENA', i_zone)
        vale_pene = mminfr(ds_contact%sdcont_defi, 'PENE_MAXI', i_zone)
        glis_maxi = mminfr(ds_contact%sdcont_defi, 'GLIS_MAXI', i_zone)
        ! l'utilisateur n'a pas renseigne glis_maxi
        if (glis_maxi .le. r8prem()) glis_maxi = 0.1*ds_contact%arete_min
        l_granglis = mminfl(ds_contact%sdcont_defi, 'GRAND_GLIS', i_zone)
!
! ----- No computation: no contact point
!
        if (l_veri) then
            goto 25
        end if
!
! ----- Loop on slave elements
!
        do i_elem_slav = 1, nb_elem_slav
!
! --------- Slave element index in contact datastructure
!
            elem_slav_indx = jdecme+i_elem_slav
!
! --------- Informations about slave element
!
            call cfnumm(ds_contact%sdcont_defi, elem_slav_indx, elem_slav_nume)
!
! --------- Number of integration points on element
!
            call mminfm(elem_slav_indx, ds_contact%sdcont_defi, 'NPTM', nb_poin_elem)
!
! --------- Get coordinates of slave element
!
            call mcomce(mesh, newgeo, elem_slav_nume, elem_slav_coor, elem_slav_type, &
                        elem_slav_nbno)
!
! --------- Get value of contact lagrangian multiplier at slave nodes
!
            call mmextm(ds_contact%sdcont_defi, cnscon, elem_slav_indx, lagr_cont_node)
!
! --------- Get value of friction lagrangian multipliers at slave nodes
!
            if (l_frot_zone) then
                call mmextm(ds_contact%sdcont_defi, cnsfr1, elem_slav_indx, lagr_fro1_node)
                if (model_ndim .eq. 3) then
                    call mmextm(ds_contact%sdcont_defi, cnsfr2, elem_slav_indx, lagr_fro2_node)
                end if
            end if
!
! --------- Loop on integration points
!
            do i_poin_elem = 1, nb_poin_elem
!
! ------------- Current master element
!
                elem_mast_nume = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+3))
!
! ------------- Get coordinates of the contact point
!
                ksipc1 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+4)
                ksipc2 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+5)
                ksipc1_old = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+19)
                ksipc2_old = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+20)
                wpg_old = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+74)
!
! ------------- Get coordinates of the projection of contact point
!
                ksipr1 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+6)
                ksipr2 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+7)
                ksipr1_old = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+22)
                ksipr2_old = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+23)
!
! ------------- Get local basis
!
                tau1(1) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+8)
                tau1(2) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+9)
                tau1(3) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+10)
                tau2(1) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+11)
                tau2(2) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+12)
                tau2(3) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+13)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+74) = &
                    v_sdcont_tabfin(ztabf*(i_cont_poin-1)+15)
!
! ------------- Store current local basis :
!               needed for previous cycling matrices and vectors computrations
!
                do hist_index = 1, 24
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+hist_index) = &
                        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+hist_index)
                end do
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+61+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+61)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+62+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+62)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+63+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+63)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+64+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+64)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+65+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+65)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+66+6) = &
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+66)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+75) = &
                    wpg_old

                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+13) = tau1(1)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+14) = tau1(2)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+15) = tau1(3)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+16) = tau2(1)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+17) = tau2(2)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+18) = tau2(3)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+19) = ksipc1
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+20) = ksipc2
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+22) = ksipr1
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+23) = ksipr2
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24) = elem_mast_nume

!
! ------------- Compute gap and contact pressure : current configuration
!
                call mmeval_prep(mesh, time_curr, model_ndim, ds_contact, &
                                 i_zone, &
                                 ksipc1, ksipc2, ksipr1, ksipr2, &
                                 tau1, tau2, &
                                 elem_slav_indx, elem_slav_nbno, &
                                 elem_slav_type, elem_slav_coor, &
                                 elem_mast_nume, &
                                 lagr_cont_node, &
                                 norm, &
                                 gap, gap_user, lagr_cont_poin, &
                                 poin_slav_coor=coor_escl_curr, poin_proj_coor=coor_proj_curr)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+61) = coor_escl_curr(1)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+62) = coor_escl_curr(2)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+63) = coor_escl_curr(3)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+64) = coor_proj_curr(1)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+65) = coor_proj_curr(2)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+66) = coor_proj_curr(3)
!
                if (l_granglis) then
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+73) = 1
                else
                    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+73) = 0
                end if

!
! ------------- Previous status and coefficients
!
                coef_cont = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2)
                coef_frot = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6)
!
! ------------- Initial bilateral contact ?
!
                l_glis_init = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+18)) .eq. 1
!
! ------------- Total gap
!
                gap = gap+gap_user
!
! ------------- Save gaps
!
                v_sdcont_jsupco(i_cont_poin) = gap_user
                v_sdcont_apjeu(i_cont_poin) = gap
!
! ------------- Excluded nodes => no contact !
!
                if (nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+19)) .eq. 1) then
                    indi_cont_curr = 0
                    goto 19
                end if
!
! ------------- Evaluate contact status
!
                call mmstac(gap, lagr_cont_poin, coef_cont, indi_cont_eval)
!
! ------------- Evaluate friction status
!
                if (l_frot_zone) then
                    call mmstaf(mesh, model_ndim, chdepd, coef_frot, &
                                elem_slav_nume, elem_slav_type, &
                                elem_slav_nbno, elem_mast_nume, ksipc1, &
                                ksipc2, ksipr1, ksipr2, lagr_fro1_node, lagr_fro2_node, &
                                tau1, tau2, norm, pres_frot, gap_user_frot, &
                                indi_frot_eval)
                else
                    indi_frot_eval = 0.d0
                    gap_user_frot(:) = 0.d0
                    pres_frot(:) = 0.d0
                end if
!
! ------------- Status treatment
!
                call mm_cycl_algo(ds_contact, l_frot_zone, &
                                  l_glis_init, type_adap, i_zone, i_cont_poin, &
                                  indi_cont_eval, indi_frot_eval, gap, lagr_cont_poin, &
                                  gap_user_frot, pres_frot, v_sdcont_cychis, v_sdcont_cyccoe, &
                                  v_sdcont_cyceta, indi_cont_curr, indi_frot_curr, loop_cont_vali, &
                                  loop_cont_conv, l_pena_frot, l_pena_cont, vale_pene, glis_maxi)
                v_sdcont_tabfin(ztabf*(i_cont_poin-1)+17) = lagr_cont_poin
!
19              continue
                if (iter_newt .ge. 2 .and. indi_cont_curr .eq. 1) then
                    do coun_bcle_geom = 1, 3
                        resi_geom = sqrt(((v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+61+6) &
                                           -coor_escl_curr(1))**2+ &
                                          (v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+62+6) &
                                           -coor_escl_curr(2))**2+ &
                                          (v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+63+6) &
                                           -coor_escl_curr(3))**2 &
                                          )+ &
                                         ((v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+64+6) &
                                           -coor_proj_curr(1))**2+ &
                                          (v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+65+6) &
                                           -coor_proj_curr(2))**2+ &
                                          (v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+66+6) &
                                           -coor_proj_curr(3))**2 &
                                          ))
                    end do
                    if (ds_contact%arete_max .gt. 1.d-15) then
                        resi_geom = resi_geom/ds_contact%arete_max
                    else
                        resi_geom = resi_geom
                    end if

                    if (resi_geom .gt. ds_contact%critere_geom .and. &
                        iter_newt .ge. 2) then
                        ds_contact%critere_geom = resi_geom
                        call codent(i_poin_elem, 'G', iptxt)
                        ds_contact%crit_geom_noeu = 'NPOINCO'//iptxt//' '
                    end if
                end if
!
! ------------- Save status
!

                v_sdcont_tabfin(ztabf*(i_cont_poin-1)+23) = indi_cont_curr
                if (l_frot_zone) then
                    v_sdcont_tabfin(ztabf*(i_cont_poin-1)+24) = indi_frot_curr
                end if
!
! ------------- Print status
!
                if (niv .ge. 2) then
                    indi_cont_prev = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+23))
                    indi_frot_prev = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+24))
                    call mmimp4(ifm, mesh, elem_slav_nume, i_poin_elem, indi_cont_prev, &
                                indi_cont_curr, indi_frot_prev, indi_frot_curr, l_frot, &
                                l_glis, gap, lagr_cont_poin)
                end if
!
! ------------- Next contact point
!
                i_cont_poin = i_cont_poin+1
                if (indi_cont_curr .eq. 1) sum_cont_press = sum_cont_press+lagr_cont_poin
            end do
        end do
25      continue
    end do
! Moyenne des pressions de contact
    sum_cont_press = sum_cont_press/nb_cont_poin
!   Critere specifique au contact en mode : RESI_CONT ne -1
    if (resi_cont .lt. r8prem()) resi_cont = 1.d-10
    resi_press_curr = abs(ds_contact%cont_pressure-abs(sum_cont_press))
    if (resi_press_curr .lt. resi_cont*abs(sum_cont_press) .and. .not. loop_cont_conv &
        .and. .not. l_loop_cont) then
        loop_cont_conv = ASTER_TRUE
    end if
    ds_contact%resi_press_glob = resi_press_curr
    ds_contact%cont_pressure = abs(sum_cont_press)
!
! - Bilateral contact management
!
    if (loop_cont_conv .and. l_exis_glis) then
        call mmglis(ds_contact)
    end if
!
! - Statistics for cycling
!
!     call mm_cycl_stat(ds_measure, ds_contact)
!
! - Propagation of coefficient
!
    l_coef_adap = ((type_adap .eq. 1) .or. (type_adap .eq. 2) .or. &
                   (type_adap .eq. 5) .or. (type_adap .eq. 6))
!     .or.  (type_adap .eq. 7) .or. (type_adap .eq. 11))
    if (l_coef_adap) then
        call mm_cycl_prop(ds_contact)
    end if
!
! - Event management for impact
!
    call mmbouc(ds_contact, 'Geom', 'Read_Counter', loop_geom_count)
    call mmbouc(ds_contact, 'Fric', 'Read_Counter', loop_fric_count)
    call mmbouc(ds_contact, 'Cont', 'Read_Counter', loop_cont_count)
    if ((iter_newt .eq. 0) .and. &
        (loop_geom_count .eq. 1) .and. (loop_fric_count .eq. 1) .and. (loop_cont_count .eq. 1)) then
        call mmeven('INI', ds_contact)
    else
        call mmeven('FIN', ds_contact)
    end if
!
! - Set loop values
!
    if (loop_cont_conv) then
        call mmbouc(ds_contact, 'Cont', 'Set_Convergence')
    else
        call mmbouc(ds_contact, 'Cont', 'Set_Divergence')
    end if
    loop_cont_vale = real(loop_cont_vali, kind=8)
    call mmbouc(ds_contact, 'Cont', 'Set_Vale', loop_vale_=loop_cont_vale)
!
! - Is contact stabilized ?
!
    if (loop_cont_vali .eq. 0) then
        ds_contact%iContStab = ds_contact%iContStab+1
    else
        ds_contact%iContStab = 0
    end if
!
! - Cleaning
!
    call jedetr(newgeo)
    call jedetr(chdepd)
    call detrsd('CHAM_NO_S', cnscon)
    call detrsd('CHAM_NO_S', cnsfr1)
    call detrsd('CHAM_NO_S', cnsfr2)
!
    call jedema()
end subroutine
