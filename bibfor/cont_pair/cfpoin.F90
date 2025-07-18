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
subroutine cfpoin(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/mminfi.h"
#include "asterfort/assert.h"
#include "asterfort/cfcorn.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmex.h"
#include "asterfort/cfnumn.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/jeexin.h"
#include "asterfort/jerazo.h"
#include "asterfort/jelira.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Fill pairing datastructure
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
    integer(kind=8) :: i_node_escl
    integer(kind=8) :: i_poin, i_node_slav, i_zone
    integer(kind=8) :: nb_poin, nb_node_slav
    integer(kind=8) :: node_slav_indx(1), node_slav_nume(1)
    integer(kind=8) :: jdecne
    real(kind=8) :: poin_coor(3)
    character(len=8) :: node_slav_name
    character(len=16) :: poin_name
    integer(kind=8) :: nb_cont_zone
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
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
!
! - Loop on contact zones
!
    i_poin = 1
    do i_zone = 1, nb_cont_zone
!
! ----- Get parameters on current zone
!
        nb_poin = mminfi(ds_contact%sdcont_defi, 'NBPT', i_zone)
        nb_node_slav = mminfi(ds_contact%sdcont_defi, 'NBNOE', i_zone)
        jdecne = mminfi(ds_contact%sdcont_defi, 'JDECNE', i_zone)
        ASSERT(nb_poin .eq. nb_node_slav)
!
! ----- Loop on contact nodes
!
        do i_node_slav = 1, nb_node_slav
!
! --------- Current contact point
!
            node_slav_indx(1) = jdecne+i_node_slav
            call cfnumn(ds_contact%sdcont_defi, 1, node_slav_indx(1), node_slav_nume(1))
!
! --------- Coordinates of contact point
!
            call cfcorn(newgeo, node_slav_nume(1), poin_coor)
            v_sdappa_poin(3*(i_poin-1)+1) = poin_coor(1)
            v_sdappa_poin(3*(i_poin-1)+2) = poin_coor(2)
            v_sdappa_poin(3*(i_poin-1)+3) = poin_coor(3)
!
! --------- Node is excluded ?
!
            call cfmmex(ds_contact%sdcont_defi, 'CONT', i_zone, node_slav_nume(1), i_node_escl)
            v_sdappa_infp(i_poin) = i_node_escl
!
! --------- Name of point
!
            node_slav_name = int_to_char8(node_slav_nume(1))
            poin_name = 'NOEUD   '//node_slav_name
            v_sdappa_noms(i_poin) = poin_name
!
! --------- Next point
!
            i_poin = i_poin+1
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
    call jerazo(sdappa_appa, 4*(i_poin-1), 1)
    call jerazo(sdappa_dist, 4*(i_poin-1), 1)
    call jerazo(sdappa_tau1, 3*(i_poin-1), 1)
    call jerazo(sdappa_tau2, 3*(i_poin-1), 1)
    call jerazo(sdappa_proj, 2*(i_poin-1), 1)
    call jerazo(sdappa_tgno, 6*(i_poin-1), 1)
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
end subroutine
