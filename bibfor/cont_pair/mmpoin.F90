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
subroutine mmpoin(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmex.h"
#include "asterfort/cfnumm.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/mcomce.h"
#include "asterfort/mmelin.h"
#include "asterfort/mmgaus.h"
#include "asterfort/mminfi.h"
#include "asterfort/mmnpoi.h"
#include "asterfort/mmnumn.h"
#include "asterfort/mmvalp.h"
#include "asterfort/jerazo.h"
#include "asterfort/jelira.h"
#include "asterfort/jeexin.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Fill pairing datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, length
    character(len=19) :: sdappa, newgeo
    character(len=24) :: sdappa_poin, sdappa_infp, sdappa_noms
    character(len=24) :: sdappa_tau1, sdappa_tau2, sdappa_proj
    character(len=24) :: sdappa_dist, sdappa_appa, sdappa_tgno, sdappa_tgel
    character(len=24) :: sdappa_mpia, sdappa_mpib, sdappa_mpic
    real(kind=8), pointer :: v_sdappa_poin(:) => null()
    integer(kind=8), pointer :: v_sdappa_infp(:) => null()
    character(len=16), pointer :: v_sdappa_noms(:) => null()
    character(len=16), pointer :: valk(:) => null()
    integer(kind=8) :: i_node_escl, i_poin, i_poin_elem, i_zone, i_elem_slav
    integer(kind=8) :: nb_elem_slav, nb_poin_elem, elem_slav_nbnode, nt_poin
    integer(kind=8) :: elem_slav_indx, elem_slav_nume, node_slav_nume
    integer(kind=8) :: jdecme
    integer(kind=8) :: type_inte
    real(kind=8) :: poin_coor(3), elem_slav_coor(27)
    real(kind=8) :: ksi1, ksi2
    character(len=8) :: elem_slav_type, elem_slav_name
    character(len=16) :: poin_name
    integer(kind=8) :: model_ndim, nb_cont_zone
!
! --------------------------------------------------------------------------------------------------
!

!
! - Pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - New geometry name
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! - Access to pairing datastructure
!
    sdappa_poin = sdappa(1:19)//'.POIN'
    sdappa_infp = sdappa(1:19)//'.INFP'
    sdappa_noms = sdappa(1:19)//'.NOMS'
    call jeveuo(sdappa_poin, 'E', vr=v_sdappa_poin)
    call jeveuo(sdappa_infp, 'E', vi=v_sdappa_infp)
    call jeveuo(sdappa_noms, 'E', vk16=v_sdappa_noms)
!
! - Get parameters
!
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
!
! - Loop on contact zones
!
    i_poin = 1
    nt_poin = 0
    do i_zone = 1, nb_cont_zone
!
! ----- Get parameters on current zone
!
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', i_zone)
        type_inte = mminfi(ds_contact%sdcont_defi, 'INTEGRATION', i_zone)
!
! ----- Loop on slave elements
!
        do i_elem_slav = 1, nb_elem_slav
!
! --------- Get current slave element
!
            elem_slav_indx = jdecme+i_elem_slav
            call cfnumm(ds_contact%sdcont_defi, elem_slav_indx, elem_slav_nume)
            elem_slav_name = int_to_char8(elem_slav_nume)
!
! --------- Get coordinates of slave element
!
            call mcomce(mesh, newgeo, elem_slav_nume, elem_slav_coor, elem_slav_type, &
                        elem_slav_nbnode)
!
! --------- Number of contact (integration) points on current slave element
!
            call mmelin(mesh, elem_slav_nume, type_inte, nb_poin_elem)
!
! --------- Loop on contact (integration) points
!
            do i_poin_elem = 1, nb_poin_elem
!
! ------------- Get current contact point
!
                call mmnumn(mesh, type_inte, elem_slav_nume, &
                            elem_slav_nbnode, i_poin_elem, node_slav_nume)
!
! ------------- Parameters of current integration point
!
                call mmgaus(elem_slav_type, type_inte, i_poin_elem, ksi1, ksi2)
!
! ------------- Coordinates of contact point
!
                call mmvalp(model_ndim, elem_slav_type, elem_slav_nbnode, 3, &
                            ksi1, ksi2, elem_slav_coor, poin_coor)
                v_sdappa_poin(3*(i_poin-1)+1) = poin_coor(1)
                v_sdappa_poin(3*(i_poin-1)+2) = poin_coor(2)
                v_sdappa_poin(3*(i_poin-1)+3) = poin_coor(3)
!
! ------------- Node is excluded ?
!
                if (node_slav_nume .gt. 0) then
                    call cfmmex(ds_contact%sdcont_defi, 'CONT', i_zone, node_slav_nume, i_node_escl)
                    v_sdappa_infp(i_poin) = i_node_escl
                end if
!
! ------------- Name of point
!
                call mmnpoi(mesh, elem_slav_name, node_slav_nume, i_poin_elem, poin_name)
                v_sdappa_noms(i_poin) = poin_name
!
! ------------- Next point
!
                i_poin = i_poin+1
                nt_poin = nt_poin+1
            end do
        end do
    end do
!
! - Pairing mpi data sutructure initialisation
!
    sdappa_appa = sdappa(1:19)//'.APPA'
    sdappa_dist = sdappa(1:19)//'.DIST'
    sdappa_tau1 = sdappa(1:19)//'.TAU1'
    sdappa_tau2 = sdappa(1:19)//'.TAU2'
    sdappa_proj = sdappa(1:19)//'.PROJ'
    sdappa_tgel = sdappa(1:19)//'.TGEL'
    sdappa_tgno = sdappa(1:19)//'.TGNO'
    call jerazo(sdappa_appa, 4*nt_poin, 1)
    call jerazo(sdappa_dist, 4*nt_poin, 1)
    call jerazo(sdappa_tau1, 3*nt_poin, 1)
    call jerazo(sdappa_tau2, 3*nt_poin, 1)
    call jerazo(sdappa_proj, 2*nt_poin, 1)
    length = cfdisi(ds_contact%sdcont_defi, 'NNOCO')
    call jerazo(sdappa_tgno, 6*length, 1)
    length = 0
    call jelira(sdappa_tgel, 'LONT', length)
    call jerazo(sdappa_tgel, length, 1)

    sdappa_mpia = sdappa(1:19)//'.MPIA'
    sdappa_mpib = sdappa(1:19)//'.MPIB'
    sdappa_mpic = sdappa(1:19)//'.MPIC'
    call jeexin(sdappa_mpia, iret)
    if (iret .eq. 0) then
        call wkvect(sdappa_mpia, 'V V K16', 1, vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call wkvect(sdappa_mpib, 'V V K16', 1, vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call wkvect(sdappa_mpic, 'V V K16', 1, vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
    else
        call jeveuo(sdappa_mpia, 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call jeveuo(sdappa_mpib, 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call jeveuo(sdappa_mpic, 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
    end if
!
    ASSERT(nt_poin .eq. cfdisi(ds_contact%sdcont_defi, 'NTPT'))
!
end subroutine
