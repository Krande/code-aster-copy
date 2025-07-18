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
subroutine mmmres(mesh, time_incr, ds_contact, disp_cumu_inst, sddisc, &
                  cnsinr, cnsper)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "jeveux.h"
#include "asterfort/apinfi.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/detrsd.h"
#include "asterfort/iseven.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcopco.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
#include "asterfort/mmmred.h"
#include "asterfort/mmmreg.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/mmfield_prep.h"
!
    character(len=8), intent(in) :: mesh
    real(kind=8), intent(in) :: time_incr
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: disp_cumu_inst
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: cnsinr
    character(len=19), intent(in) :: cnsper
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Post-treatment
!
! Continue method - Prepare post-treatment fields
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  time_incr        : time increment
! In  ds_contact       : datastructure for contact management
! In  disp_cumu_inst   : displacement increment from beginning of current time
! In  sddisc           : datastructure for discretization
! In  cnsinr           : nodal field (CHAM_NO_S) for CONT_NOEU
! In  cnsper           : nodal field (CHAM_NO_S) to save percussions
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: eps = 1.d-6
    integer(kind=8) :: i_cont_poin, i_zone, i_elem_slav, i_poin_elem
    integer(kind=8) :: nb_dof, nb_poin_elem, nb_elem_slav, nb_cont_zone, model_ndim, nt_cont_poin
    integer(kind=8) :: ztabf, zresu, zperc
    integer(kind=8) :: node_slav_nume, elem_mast_nume, elem_slav_indx
    integer(kind=8) :: jdecme, indi_cont, pair_type
    real(kind=8) :: gli, gli1, gli2
    real(kind=8) :: rn, rnx, rny, rnz
    real(kind=8) :: rtax, rtay, rtaz
    real(kind=8) :: rtgx, rtgy, rtgz
    real(kind=8) :: r, rx, ry, rz
    real(kind=8) :: imp, impx, impy, impz
    real(kind=8) :: node_status, lagsf, lagf(2)
    real(kind=8) :: ksipr1, ksipr2, proj(3)
    character(len=19) :: disp_cumu_s
    character(len=19) :: cneltc_s, cneltf_s
    character(len=19) :: cneltc, cneltf
    character(len=19) :: newgeo, sdappa
    character(len=24) :: sdcont_tabfin, sdcont_apjeu
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    real(kind=8), pointer :: v_sdcont_apjeu(:) => null()
    real(kind=8) :: valras
    aster_logical :: l_frot_zone, l_veri, lcolli, laffle, l_frot
    real(kind=8), pointer :: v_cnsper_cnsv(:) => null()
    real(kind=8), pointer :: v_cnsinr_cnsv(:) => null()
    aster_logical, pointer :: v_cnsper_cnsl(:) => null()
    aster_logical, pointer :: v_cnsinr_cnsl(:) => null()
    real(kind=8), pointer :: v_cneltc(:) => null()
    real(kind=8), pointer :: v_cneltf(:) => null()
    real(kind=8), pointer :: v_slav_slide(:) => null()
    real(kind=8), pointer :: v_mast_slide(:) => null()
    real(kind=8), pointer :: v_disp_cumu(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    cneltc = '&&MMMRES.CONT'
    cneltc_s = '&&MMMRES.FCTCN'
    cneltf = '&&MMMRES.FROT'
    cneltf_s = '&&MMMRES.FROTCN'
    disp_cumu_s = '&&MMMRES.DEPCN'
    lagf = 0.d0
!
! - Get parameters
!
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROTTEMENT')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    nt_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
!
! - Collision
!
    laffle = ASTER_FALSE
    valras = 1.d-3
    call iseven(sddisc, FAIL_EVT_COLLISION, lcolli)
!
! - Acces to contact objects
!
    ztabf = cfmmvd('ZTABF')
    zperc = cfmmvd('ZPERC')
    zresu = cfmmvd('ZRESU')
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    sdcont_apjeu = ds_contact%sdcont_solv(1:14)//'.APJEU'
    call jeveuo(sdcont_tabfin, 'E', vr=v_sdcont_tabfin)
    call jeveuo(sdcont_apjeu, 'E', vr=v_sdcont_apjeu)
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - Geometric update
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! - Get fields
!
    cneltc = ds_contact%cneltc
    cneltf = ds_contact%cneltf
!
! - Prepare nodal fields (select right components)
!
    call mmmred(model_ndim, l_frot, disp_cumu_inst, disp_cumu_s, nb_dof)
    call mmfield_prep(cneltc, cneltc_s, &
                      l_sort_=.true._1, nb_cmp_=model_ndim, &
                      list_cmp_=['DX      ', 'DY      ', 'DZ      '])
    if (l_frot) then
        call mmfield_prep(cneltf, cneltf_s, &
                          l_sort_=.true._1, nb_cmp_=model_ndim, &
                          list_cmp_=['DX      ', 'DY      ', 'DZ      '])
    end if
!
! - Access to fields
!
    call jeveuo(disp_cumu_s(1:19)//'.CNSV', 'L', vr=v_disp_cumu)
    call jeveuo(cnsinr(1:19)//'.CNSV', 'E', vr=v_cnsinr_cnsv)
    call jeveuo(cnsinr(1:19)//'.CNSL', 'E', vl=v_cnsinr_cnsl)
    call jeveuo(cnsper(1:19)//'.CNSV', 'E', vr=v_cnsper_cnsv)
    call jeveuo(cnsper(1:19)//'.CNSL', 'E', vl=v_cnsper_cnsl)
    call jeveuo(cneltc_s(1:19)//'.CNSV', 'L', vr=v_cneltc)
    if (l_frot) then
        call jeveuo(cneltf_s(1:19)//'.CNSV', 'L', vr=v_cneltf)
    end if
!
! - Compute slide
!
    AS_ALLOCATE(vr=v_slav_slide, size=2*nt_cont_poin)
    AS_ALLOCATE(vr=v_mast_slide, size=2*nt_cont_poin)
    call mmmreg(mesh, ds_contact, v_disp_cumu, nb_dof, &
                v_slav_slide, v_mast_slide)
!
! - Loop on contact zones
!
    i_cont_poin = 1
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters of zone
!
        l_veri = mminfl(ds_contact%sdcont_defi, 'VERIF', i_zone)
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', i_zone)
        l_frot_zone = mminfl(ds_contact%sdcont_defi, 'FROTTEMENT_ZONE', i_zone)
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
! --------- Number of integration points on element
!
            call mminfm(elem_slav_indx, ds_contact%sdcont_defi, 'NPTM', nb_poin_elem)
!
! --------- Loop on integration points
!
            do i_poin_elem = 1, nb_poin_elem
!
                gli = 0.d0
                gli1 = 0.d0
                gli2 = 0.d0
                rtax = 0.d0
                rtay = 0.d0
                rtaz = 0.d0
                rtgx = 0.d0
                rtgy = 0.d0
                rtgz = 0.d0
                rn = 0.d0
                rnx = 0.d0
                rny = 0.d0
                rnz = 0.d0
                node_status = 0.d0
!
! ------------- Get pairing
!
                call apinfi(sdappa, 'APPARI_TYPE', i_poin_elem, pair_type)
                if (pair_type .lt. 0) then
                    node_status = -1.d0
                end if
!
! ------------- Get slave node index
!
                node_slav_nume = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+25))
                if (node_slav_nume .le. 0) then
                    goto 99
                end if
!
! ------------- Contact status
!
                indi_cont = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+23))
!
! ------------- Projection
!
                ksipr1 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+6)
                ksipr2 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+7)
                elem_mast_nume = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+3))
                call mcopco(mesh, newgeo, model_ndim, elem_mast_nume, ksipr1, &
                            ksipr2, proj)
!
! ------------- Get information for contact
!
                if (indi_cont .eq. 1) then
!
                    node_status = 1.d0
!
! ------------- Get nodal reaction for contact
!
                    if (model_ndim .eq. 3) then
                        rnx = v_cneltc(3*(node_slav_nume-1)+1)
                        rny = v_cneltc(3*(node_slav_nume-1)+2)
                        rnz = v_cneltc(3*(node_slav_nume-1)+3)
                        if (abs(rnx) > 10.d150 .or. abs(rny) > 10.d150 .or. abs(rnz) > 10.d150) then
                            rn = 10.d150
                        else
                            rn = sqrt(rnx**2+rny**2+rnz**2)
                        end if
                    else if (model_ndim .eq. 2) then
                        rnx = v_cneltc(2*(node_slav_nume-1)+1)
                        rny = v_cneltc(2*(node_slav_nume-1)+2)
                        if (abs(rnx) > 10.d150 .or. abs(rny) > 10.d150) then
                            rn = 10.d150
                        else
                            rn = sqrt(rnx**2+rny**2)
                        end if
                    else
                        ASSERT(ASTER_FALSE)
                    end if
!
! ----------------- Very near contact
!
                    if (rn .le. valras) then
                        laffle = ASTER_TRUE
                    end if
!
! ----------------- Friction
!
                    if (l_frot_zone) then
!
! --------------------- Compute slides
!
                        if (model_ndim .eq. 3) then
                            gli1 = v_slav_slide(2*(i_cont_poin-1)+1)- &
                                   v_mast_slide(2*(i_cont_poin-1)+1)
                            gli2 = v_slav_slide(2*(i_cont_poin-1)+2)- &
                                   v_mast_slide(2*(i_cont_poin-1)+2)
                            if (abs(gli1) > 10.d150 .or. abs(gli2) > 10.d150) then
                                gli = 10.d150
                            else
                                gli = sqrt(gli1**2+gli2**2)
                            end if
                        else if (model_ndim .eq. 2) then
                            gli1 = v_slav_slide(i_cont_poin)- &
                                   v_mast_slide(i_cont_poin)
                            gli = abs(gli1)
                        else
                            ASSERT(ASTER_FALSE)
                        end if
!
! --------------------- Friction Lagrange
!
                        if (model_ndim .eq. 3) then
!
! --------------------- Test to prevent FPE
!
                            lagf(1) = v_disp_cumu(nb_dof*(node_slav_nume-1)+5)
                            lagf(2) = v_disp_cumu(nb_dof*(node_slav_nume-1)+6)
                            if (maxval(abs(lagf)) > 10.d50) then
                                lagsf = 10.d50
                            else
                                lagsf = norm2(lagf)
                            end if
                        else if (model_ndim .eq. 2) then
                            lagsf = abs(v_disp_cumu(nb_dof*(node_slav_nume-1)+4))
                        else
                            ASSERT(ASTER_FALSE)
                        end if
!
! --------------------- Stick or slide ?
!
                        if (lagsf .ge. 0.999d0) then
                            node_status = 2.d0
                            if (model_ndim .eq. 3) then
                                rtgx = v_cneltf(3*(node_slav_nume-1)+1)
                                rtgy = v_cneltf(3*(node_slav_nume-1)+2)
                                rtgz = v_cneltf(3*(node_slav_nume-1)+3)
                            else if (model_ndim .eq. 2) then
                                rtgx = v_cneltf(2*(node_slav_nume-1)+1)
                                rtgy = v_cneltf(2*(node_slav_nume-1)+2)
                            else
                                ASSERT(ASTER_FALSE)
                            end if
                        else
                            node_status = 1.d0
                            if (model_ndim .eq. 3) then
                                rtax = v_cneltf(3*(node_slav_nume-1)+1)
                                rtay = v_cneltf(3*(node_slav_nume-1)+2)
                                rtaz = v_cneltf(3*(node_slav_nume-1)+3)
                            else if (model_ndim .eq. 2) then
                                rtax = v_cneltf(2*(node_slav_nume-1)+1)
                                rtay = v_cneltf(2*(node_slav_nume-1)+2)
                            else
                                ASSERT(ASTER_FALSE)
                            end if
                        end if
                    else
                        lagsf = 0.d0
                        node_status = 2.d0
                    end if
                end if
!
! ------------- Total reaction
!
                rx = rnx+rtax+rtgx
                ry = rny+rtay+rtgy
                rz = rnz+rtaz+rtgz
! ------------- To prevent FPE
                if (abs(rx) > 10.d150 .or. abs(ry) > 10.d150 .or. abs(rz) > 10.d150) then
                    r = 10.d150
                else
                    r = sqrt(rx**2+ry**2+rz**2)
                end if
!
! ------------- Percussion
!
                if (r .le. eps) then
                    imp = 0.d0
                    impx = 0.d0
                    impy = 0.d0
                    impz = 0.d0
                else
                    imp = v_cnsper_cnsv(zperc*(node_slav_nume-1)+1)+r*time_incr
                    impx = v_cnsper_cnsv(zperc*(node_slav_nume-1)+2)+rx*time_incr
                    impy = v_cnsper_cnsv(zperc*(node_slav_nume-1)+3)+ry*time_incr
                    impz = v_cnsper_cnsv(zperc*(node_slav_nume-1)+4)+rz*time_incr
                end if
!
! ------------- Save in CONT_NOEU field
!
                v_cnsinr_cnsv(zresu*(node_slav_nume-1)+1) = node_status
                v_cnsinr_cnsv(zresu*(node_slav_nume-1)+2) = -v_sdcont_apjeu(i_cont_poin)
                v_cnsinr_cnsl(zresu*(node_slav_nume-1)+1) = ASTER_TRUE
                v_cnsinr_cnsl(zresu*(node_slav_nume-1)+2) = ASTER_TRUE
                if (model_ndim .eq. 3) then
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+3) = rn
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+4) = rnx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+5) = rny
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+6) = rnz
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+7) = gli1
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+8) = gli2
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+9) = gli
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+10) = rtax
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+11) = rtay
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+12) = rtaz
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+13) = rtgx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+14) = rtgy
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+15) = rtgz
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+16) = rx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+17) = ry
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+18) = rz
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+19) = r
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+21) = imp
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+22) = impx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+23) = impy
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+24) = impz
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+28) = proj(1)
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+29) = proj(2)
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+30) = proj(3)
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+3) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+4) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+5) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+6) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+7) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+8) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+9) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+10) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+11) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+12) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+13) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+14) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+15) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+16) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+17) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+18) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+19) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+21) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+22) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+23) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+24) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+28) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+29) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+30) = ASTER_TRUE
                else if (model_ndim .eq. 2) then
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+3) = rn
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+4) = rnx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+5) = rny
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+7) = gli1
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+9) = gli
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+10) = rtax
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+11) = rtay
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+13) = rtgx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+14) = rtgy
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+16) = rx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+17) = ry
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+19) = r
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+21) = imp
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+22) = impx
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+23) = impy
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+28) = proj(1)
                    v_cnsinr_cnsv(zresu*(node_slav_nume-1)+29) = proj(2)
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+3) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+4) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+5) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+7) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+9) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+10) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+11) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+13) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+14) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+16) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+17) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+19) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+21) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+22) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+23) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+28) = ASTER_TRUE
                    v_cnsinr_cnsl(zresu*(node_slav_nume-1)+29) = ASTER_TRUE
                else
                    ASSERT(ASTER_FALSE)
                end if
!
! ------------- Save in percussion field
!
                v_cnsper_cnsv(zperc*(node_slav_nume-1)+1) = imp
                v_cnsper_cnsv(zperc*(node_slav_nume-1)+2) = impx
                v_cnsper_cnsv(zperc*(node_slav_nume-1)+3) = impy
                v_cnsper_cnsl(zperc*(node_slav_nume-1)+1) = ASTER_TRUE
                v_cnsper_cnsl(zperc*(node_slav_nume-1)+2) = ASTER_TRUE
                v_cnsper_cnsl(zperc*(node_slav_nume-1)+3) = ASTER_TRUE
!
                if (model_ndim .eq. 3) then
                    v_cnsper_cnsv(zperc*(node_slav_nume-1)+4) = impz
                    v_cnsper_cnsl(zperc*(node_slav_nume-1)+4) = ASTER_TRUE
                end if
99              continue
!
! ------------- Next contact point
!
                i_cont_poin = i_cont_poin+1
            end do
        end do
25      continue
    end do
!
! - Clean
!
    call jedetr(cneltc)
    call detrsd('CHAMP', cneltc_s)
    call jedetr(cneltf)
    call detrsd('CHAMP', cneltf_s)
    call detrsd('CHAMP', disp_cumu_s)
    AS_DEALLOCATE(vr=v_slav_slide)
    AS_DEALLOCATE(vr=v_mast_slide)
!
! - Alarm for COLLISION
!
    if (laffle .and. lcolli) then
        call utmess('A', 'CONTACT3_98')
    end if
!
! - Alarm for large penetration criteria
!
    if (nint(ds_contact%continue_pene) .eq. 1) then
        call utmess('A', 'CONTACT3_99', nr=2, valr=[ds_contact%arete_min, ds_contact%arete_max])
    end if
!
! - Alarm for penetration criteria
!
!     if (ds_contact%calculated_penetration .le. 1.d-99) then
!         call utmess('A', 'CONTACT3_97')
!     endif
!
    call jedema()
end subroutine
