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

subroutine apimpr_l(ifm, mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cfdisi.h"
#include "asterfort/codent.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/mminfi.h"
#include "asterfort/get_patch_info.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: ifm
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Segment to segment - Debug print
!
! --------------------------------------------------------------------------------------------------
!
! In  ifm              : unit for message
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: knuzo
    integer(kind=8) :: nb_elem_mast, nb_elem_slav, nb_node_slav, nb_node_mast, nb_cont_pair
    integer(kind=8) :: nb_elem_patch
    integer(kind=8) :: list_elem(5)
    character(len=24) :: sdappa_mast, sdappa_slav
    integer(kind=8), pointer :: v_sdappa_mast(:) => null()
    integer(kind=8), pointer :: v_sdappa_slav(:) => null()
    integer(kind=8), pointer :: v_mesh_comapa(:) => null()
    integer(kind=8), pointer :: v_mesh_lpatch(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    integer(kind=8) :: i_cont_zone, i_patch, i_node, i_elem_slav, i_cont_pair, i_elem_patch, i_elem
    integer(kind=8) :: patch_jdec, patch_indx, patch_nume
    integer(kind=8) :: nb_cont_zone, elem_nbnode
    integer(kind=8) :: nt_patch, nt_node_slav, nt_node_mast
    integer(kind=8) :: node_nume(9), elem_slav_nume, elem_mast_nume, elem_nume
    character(len=8) :: node_name(9), elem_name(9), elem_slav_name, elem_mast_name
    real(kind=8) :: coefint, gap, weight, weight_total
    character(len=19) :: sdappa
    character(len=24) :: sdappa_gapi, sdappa_coef, sdappa_poid
    real(kind=8), pointer :: v_sdappa_gapi(:) => null()
    real(kind=8), pointer :: v_sdappa_coef(:) => null()
    real(kind=8), pointer :: v_sdappa_poid(:) => null()
    character(len=24) :: sdappa_apli
    integer(kind=8), pointer :: v_sdappa_apli(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    weight_total = 0.d0
!
! - Get parameters
!
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    nt_node_slav = cfdisi(ds_contact%sdcont_defi, 'NTNOE')
    nt_node_mast = cfdisi(ds_contact%sdcont_defi, 'NTNOM')
    nb_cont_pair = ds_contact%nb_cont_pair
    nt_patch = ds_contact%nt_patch
!
! - Access to mesh (patches)
!
    call jeveuo(jexnum(mesh//'.PATCH', 1), 'L', vi=v_mesh_lpatch)
    call jeveuo(mesh//'.COMAPA', 'L', vi=v_mesh_comapa)
!
! - Access to pairing datastructures
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
    sdappa_gapi = sdappa(1:19)//'.GAPI'
    sdappa_coef = sdappa(1:19)//'.COEF'
    sdappa_poid = sdappa(1:19)//'.POID'
    sdappa_apli = sdappa(1:19)//'.APLI'
    call jeveuo(sdappa_gapi, 'L', vr=v_sdappa_gapi)
    call jeveuo(sdappa_coef, 'L', vr=v_sdappa_coef)
    call jeveuo(sdappa_poid, 'L', vr=v_sdappa_poid)
    call jeveuo(sdappa_apli, 'L', vi=v_sdappa_apli)
!
! - Patchs
!
    write (ifm, 10)
    write (ifm, 101) nt_patch
    do i_patch = 1, nt_patch
        patch_indx = i_patch
        call get_patch_info(sdappa, patch_indx, nb_elem_patch, list_elem)
        do i_elem_patch = 1, nb_elem_patch
            elem_nume = list_elem(i_elem_patch)
            elem_name(i_elem_patch) = int_to_char8(elem_nume)
        end do
        write (ifm, 302) i_patch, (elem_name(i_elem), i_elem=1, nb_elem_patch)
        coefint = v_sdappa_coef(patch_indx)
        gap = v_sdappa_gapi(patch_indx)
        weight = v_sdappa_poid(patch_indx)
        if (isnan(gap)) then
            write (ifm, 303) coefint, weight
        else
            write (ifm, 304) gap, coefint, weight
        end if
    end do
302 format(' <Pairing> . Patch  : ', i5, ' with elements ', 5(a8, 1x))
303 format(' <Pairing> .. No Gap : Inter. coef: ', e12.5, &
           ' - Inter weight:', e12.5)
304 format(' <Pairing> .. Gap    :', e12.5, ' - Inter. coef: ', e12.5, &
           ' - Inter weight:', e12.5)
!
! - Contact zones
!
    write (ifm, 10)
    write (ifm, 100) nb_cont_zone
    do i_cont_zone = 1, nb_cont_zone
!
! ----- Get parameters of current zone
!
        ASSERT(i_cont_zone .le. 9)
        patch_jdec = v_mesh_lpatch(2*(i_cont_zone-1)+1)-1
        nb_node_mast = mminfi(ds_contact%sdcont_defi, 'NBNOM', i_cont_zone)
        nb_node_slav = mminfi(ds_contact%sdcont_defi, 'NBNOE', i_cont_zone)
        nb_elem_mast = mminfi(ds_contact%sdcont_defi, 'NBMAM', i_cont_zone)
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_cont_zone)
        write (ifm, 102) i_cont_zone
        write (ifm, 103) nb_elem_slav
        write (ifm, 104) nb_elem_mast
        write (ifm, 105) nb_node_slav
        write (ifm, 106) nb_node_mast
!
! ----- Get master and slaves elements of current zone
!
        call codent(i_cont_zone-1, 'G', knuzo)
        sdappa_mast = sdappa(1:19)//'.MS'//knuzo(1:2)
        sdappa_slav = sdappa(1:19)//'.EC'//knuzo(1:2)
        call jeveuo(sdappa_mast, 'L', vi=v_sdappa_mast)
        call jeveuo(sdappa_slav, 'L', vi=v_sdappa_slav)
!
! ----- Loop on slave elements
!
        do i_elem_slav = 1, nb_elem_slav
!
! --------- Get current slave element
!
            elem_slav_nume = v_sdappa_slav(i_elem_slav)
            elem_slav_name = int_to_char8(elem_slav_nume)
            call jeveuo(jexnum(mesh//'.CONNEX', elem_slav_nume), 'L', vi=v_mesh_connex)
            call jelira(jexnum(mesh//'.CONNEX', elem_slav_nume), 'LONMAX', elem_nbnode)
            node_name(:) = ' '
            do i_node = 1, elem_nbnode
                node_nume(i_node) = v_mesh_connex(i_node)
                node_name(i_node) = int_to_char8(node_nume(i_node))
            end do
            write (ifm, 109) elem_slav_name, (node_name(i_node), i_node=1, elem_nbnode)
!
! --------- Get current patch
!
            patch_indx = v_mesh_comapa(elem_slav_nume)
            patch_nume = patch_indx+1-patch_jdec
!
! --------- Get parameters
!
            coefint = v_sdappa_coef(patch_nume)
            gap = v_sdappa_gapi(patch_nume)
            weight = v_sdappa_poid(patch_nume)

            if (isnan(gap)) then
                write (ifm, 120)
            else
                write (ifm, 121) gap, coefint, weight
                weight_total = weight_total+weight
            end if
        end do
    end do
!
! - Contact pairs
!
    write (ifm, 10)
    write (ifm, 200) nb_cont_pair
    if (nb_cont_pair .ne. 0) then
        call jeveuo(sdappa_apli, 'L', vi=v_sdappa_apli)
    end if
    do i_cont_pair = 1, nb_cont_pair
        elem_slav_nume = v_sdappa_apli(3*(i_cont_pair-1)+1)
        elem_mast_nume = v_sdappa_apli(3*(i_cont_pair-1)+2)
        elem_slav_name = int_to_char8(elem_slav_nume)
        elem_mast_name = int_to_char8(elem_mast_nume)
        write (ifm, 201) i_cont_pair, elem_slav_name, elem_mast_name
    end do
    write (ifm, 202) weight_total
!
10  format(' <Pairing> **************************************************************************')
100 format(' <Pairing> Number of zones         : ', i6)
101 format(' <Pairing> Total number of patches : ', i6)
102 format(' <Pairing> . Zone : ', i6)
103 format(' <Pairing> .. Total number of slave elements  : ', i6)
104 format(' <Pairing> .. Total number of master elements : ', i6)
105 format(' <Pairing> .. Total number of slave nodes     : ', i6)
106 format(' <Pairing> .. Total number of master nodes    : ', i6)
109 format(' <Pairing> ... Slave element : ', a8, ' with nodes ', 9(a8, 1x))
120 format(' <Pairing> ... Not paired')
121 format(' <Pairing> ... Paired with gap :', e12.5, ' - Inter. coef: ', e12.5, &
           ' - Inter weight:', e12.5)
200 format(' <Pairing> Number of contact pairs : ', i6)
201 format(' <Pairing> . Contact pair : ', i6, ' - Slave: ', a8, ' - Master:', a8)
202 format(' <Pairing> . Total surface of patchs: ', e12.5)
!
end subroutine
