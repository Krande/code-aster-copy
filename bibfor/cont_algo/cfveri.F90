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

subroutine cfveri(mesh, ds_contact, time_curr, nt_ncomp_poin, &
                  v_ncomp_jeux, v_ncomp_loca, v_ncomp_enti, v_ncomp_zone)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/apcopt.h"
#include "asterfort/apinfi.h"
#include "asterfort/apinfr.h"
#include "asterfort/apvect.h"
#include "asterfort/assert.h"
#include "asterfort/cfcoor.h"
#include "asterfort/cfcorn.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdist.h"
#include "asterfort/cfnewj.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfnumn.h"
#include "asterfort/cftanr.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mmnorm.h"
#include "asterfort/mmnpoi.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8), intent(in) :: time_curr
    integer(kind=8), intent(in) :: nt_ncomp_poin
    real(kind=8), pointer :: v_ncomp_jeux(:)
    integer(kind=8), pointer :: v_ncomp_loca(:)
    character(len=16), pointer :: v_ncomp_enti(:)
    integer(kind=8), pointer :: v_ncomp_zone(:)

!
! -----------------------------------------------------------------------------------------------
!
! Contact - Post-treatment for no computation methods
!
! Method discrete - Evaluate
!
! -----------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  time_curr        : current time
! In  nt_ncomp_poin    : number of points in no-computation mode
! In  v_ncomp_jeux     : pointer to save gaps
! In  v_ncomp_loca     : pointer to save index of node
! In  v_ncomp_enti     : pointer to save name of entities
! In  v_ncomp_zone     : pointer to save contact zone index
!
! -----------------------------------------------------------------------------------------------
!
    character(len=19) :: newgeo, sdappa
    integer(kind=8) :: pair_type, pair_enti
    integer(kind=8) :: jdecne
    integer(kind=8) :: posmae, elem_mast_nume, node_slav_indx(1), elem_mast_indx, node_slav_nume(1)
    integer(kind=8) :: i_zone, i_poin, i_cont_node, i_ncomp_poin
    integer(kind=8) :: model_ndim, nb_cont_zone
    integer(kind=8) :: nb_poin, nb_cont_node, nt_ncomp_poin0
    real(kind=8) :: node_coor_proj(3), poin_coor(3)
    real(kind=8) :: tau1m(3), tau2m(3)
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: norm(3), noor
    real(kind=8) :: ksipr1, ksipr2
    real(kind=8) :: r8bid
    real(kind=8) :: gap, gap_user
    character(len=8) :: node_slav_name, elem_mast_name, k8bla
    character(len=16) :: poin_name, enti_name
    aster_logical :: l_veri
!
! -----------------------------------------------------------------------------------------------
!
    node_slav_indx = 0
    i_ncomp_poin = 1
    k8bla = ' '
!
! - Pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - New geometry name
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! - Parameters
!
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
!
! - Loop on contact zones
!
    i_poin = 1
    nt_ncomp_poin0 = 0
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters of zone
!
        nb_cont_node = mminfi(ds_contact%sdcont_defi, 'NBPT', i_zone)
        jdecne = mminfi(ds_contact%sdcont_defi, 'JDECNE', i_zone)
        l_veri = mminfl(ds_contact%sdcont_defi, 'VERIF', i_zone)
!
! ----- Computation: no evaluate (see cfresu)
!
        if (.not. l_veri) then
            nb_poin = mminfi(ds_contact%sdcont_defi, 'NBPC', i_zone)
            i_poin = i_poin+nb_poin
            goto 25
        end if
!
! ----- Loop on nodes
!
        do i_cont_node = 1, nb_cont_node
!
! --------- Current slave node
!
            node_slav_indx(1) = jdecne+i_cont_node
            call cfnumn(ds_contact%sdcont_defi, 1, node_slav_indx(1), node_slav_nume(1))
            node_slav_name = int_to_char8(node_slav_nume(1))
!
! --------- Parameters from pairing
!
            call apinfi(sdappa, 'APPARI_TYPE', i_poin, pair_type)
            call apinfi(sdappa, 'APPARI_ENTITE', i_poin, pair_enti)
            call apinfr(sdappa, 'APPARI_PROJ_KSI1', i_poin, ksipr1)
            call apinfr(sdappa, 'APPARI_PROJ_KSI2', i_poin, ksipr2)
            call apvect(sdappa, 'APPARI_TAU1', i_poin, tau1m)
            call apvect(sdappa, 'APPARI_TAU2', i_poin, tau2m)
!
! --------- Coordinates of point
!
            call apcopt(sdappa, i_poin, poin_coor)
!
! --------- Name of point
!
            call mmnpoi(mesh, k8bla, node_slav_nume(1), i_cont_node, poin_name)
!
! --------- Depending on pairing type
!
            if (pair_type .eq. 2) then
!
! ------------- Master element
!
                elem_mast_indx = pair_enti
                call cfnumm(ds_contact%sdcont_defi, elem_mast_indx, elem_mast_nume)
                elem_mast_name = int_to_char8(elem_mast_nume)
                enti_name = elem_mast_name
!
! ------------- Coordinates of projection
!
                call cfcoor(mesh, ds_contact%sdcont_defi, newgeo, elem_mast_indx, ksipr1, &
                            ksipr2, node_coor_proj)
!
! ------------- Define new local basis
!
                call cftanr(mesh, model_ndim, ds_contact, i_zone, &
                            node_slav_indx(1), 'MAIL', elem_mast_indx, elem_mast_nume, ksipr1, &
                            ksipr2, tau1m, tau2m, tau1, tau2)
                call mmnorm(model_ndim, tau1, tau2, norm, noor)
                if (noor .le. r8prem()) then
                    call utmess('F', 'CONTACT3_26', sk=node_slav_name)
                end if
!
! ------------- Compute gap
!
                call cfnewj(model_ndim, poin_coor, node_coor_proj, norm, gap)
                call cfdist(ds_contact, i_zone, posmae, poin_coor, time_curr, &
                            gap_user, node_slav_indx_=node_slav_indx(1))
                gap = gap+gap_user
            else if (pair_type .eq. 1) then
!
! ------------- Master node
!
                node_slav_indx(1) = pair_enti
                call cfnumn(ds_contact%sdcont_defi, 1, node_slav_indx(1), node_slav_nume(1))
                node_slav_name = int_to_char8(node_slav_nume(1))
                enti_name = node_slav_name
!
! ------------- Coordinate of master node
!
                call cfcorn(newgeo, node_slav_nume(1), node_coor_proj)
!
! ------------- Define new local basis
!
                call cftanr(mesh, model_ndim, ds_contact, i_zone, &
                            node_slav_indx(1), 'NOEU', node_slav_indx(1), node_slav_nume(1), &
                            r8bid, r8bid, tau1m, tau2m, tau1, tau2)
                call mmnorm(model_ndim, tau1, tau2, norm, noor)
                if (noor .le. r8prem()) then
                    call utmess('F', 'CONTACT3_26', sk=node_slav_name)
                end if
!
! ------------- Compute gap
!
                call cfnewj(model_ndim, poin_coor, node_coor_proj, norm, gap)
                call cfdist(ds_contact, i_zone, posmae, poin_coor, time_curr, &
                            gap_user, node_slav_indx_=node_slav_indx(1))
                gap = gap+gap_user
!
            else if (pair_type .eq. -1) then
                enti_name = 'EXCLU'
                gap = r8vide()
            else if (pair_type .eq. -2) then
                enti_name = 'EXCLU'
                gap = r8vide()
            else if (pair_type .eq. -3) then
                enti_name = 'EXCLU'
                gap = r8vide()
            else
                ASSERT(.false.)
            end if
!
! --------- Save
!
            v_ncomp_jeux(i_ncomp_poin) = gap
            v_ncomp_loca(i_ncomp_poin) = node_slav_nume(1)
            v_ncomp_zone(i_ncomp_poin) = i_zone
            v_ncomp_enti(2*(i_ncomp_poin-1)+1) = poin_name
            v_ncomp_enti(2*(i_ncomp_poin-1)+2) = enti_name
!
! --------- Next points
!
            i_ncomp_poin = i_ncomp_poin+1
            nt_ncomp_poin0 = nt_ncomp_poin0+1
            i_poin = i_poin+1
!
        end do
25      continue
    end do
!
    ASSERT(nt_ncomp_poin0 .eq. nt_ncomp_poin)
!
end subroutine
