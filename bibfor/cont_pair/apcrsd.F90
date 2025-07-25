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
subroutine apcrsd(ds_contact, sdappa, &
                  nt_poin, nb_cont_elem, nb_cont_node, &
                  nt_elem_node, nb_node_mesh)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/cfnben.h"
#include "asterfort/jecroc.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(in) :: nt_poin
    integer(kind=8), intent(in) :: nb_cont_elem
    integer(kind=8), intent(in) :: nb_cont_node
    integer(kind=8), intent(in) :: nt_elem_node
    integer(kind=8), intent(in) :: nb_node_mesh
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Create datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  sdappa           : name of pairing datastructure
! In  nt_poin          : total number of points (contact and non-contact)
! In  nb_cont_elem     : total number of contact elements
! In  nb_cont_node     : total number of contact nodes
! In  nt_elem_node     : total number of nodes at all contact elements
! In  nb_node_mesh     : number of nodes in mesh
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_cont_elem, longt, elem_indx, longc, elem_nbnode
    character(len=24) :: sdappa_poin
    real(kind=8), pointer :: v_sdappa_poin(:) => null()
    character(len=24) :: sdappa_infp
    integer(kind=8), pointer :: v_sdappa_infp(:) => null()
    character(len=24) :: sdappa_noms
    character(len=16), pointer :: v_sdappa_noms(:) => null()
    character(len=24) :: sdappa_appa
    integer(kind=8), pointer :: v_sdappa_appa(:) => null()
    character(len=24) :: sdappa_dist
    real(kind=8), pointer :: v_sdappa_dist(:) => null()
    character(len=24) :: sdappa_tau1
    real(kind=8), pointer :: v_sdappa_tau1(:) => null()
    character(len=24) :: sdappa_tau2
    real(kind=8), pointer :: v_sdappa_tau2(:) => null()
    character(len=24) :: sdappa_proj
    real(kind=8), pointer :: v_sdappa_proj(:) => null()
    character(len=24) :: sdappa_tgno, sdappa_tgel
    real(kind=8), pointer :: v_sdappa_tgno(:) => null()
    character(len=24) :: sdappa_verk
    character(len=8), pointer :: v_sdappa_verk(:) => null()
    character(len=24) :: sdappa_vera
    real(kind=8), pointer :: v_sdappa_vera(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Datastructure for pairing results
!
    sdappa_appa = sdappa(1:19)//'.APPA'
    call wkvect(sdappa_appa, 'V V I', 4*nt_poin, vi=v_sdappa_appa)
!
! - Datastructure for distances and local basis
!
    sdappa_dist = sdappa(1:19)//'.DIST'
    sdappa_tau1 = sdappa(1:19)//'.TAU1'
    sdappa_tau2 = sdappa(1:19)//'.TAU2'
    call wkvect(sdappa_dist, 'V V R', 4*nt_poin, vr=v_sdappa_dist)
    call wkvect(sdappa_tau1, 'V V R', 3*nt_poin, vr=v_sdappa_tau1)
    call wkvect(sdappa_tau2, 'V V R', 3*nt_poin, vr=v_sdappa_tau2)
!
! - Datastructure for projection points
!
    sdappa_proj = sdappa(1:19)//'.PROJ'
    call wkvect(sdappa_proj, 'V V R', 2*nt_poin, vr=v_sdappa_proj)
!
! - Datastructure for coordinates of points
!
    sdappa_poin = sdappa(1:19)//'.POIN'
    call wkvect(sdappa_poin, 'V V R', 3*nt_poin, vr=v_sdappa_poin)
!
! - Datastructure for informations about points
!
    sdappa_infp = sdappa(1:19)//'.INFP'
    call wkvect(sdappa_infp, 'V V I', nt_poin, vi=v_sdappa_infp)
!
! - Datastructure for name of contact points
!
    sdappa_noms = sdappa(1:19)//'.NOMS'
    call wkvect(sdappa_noms, 'V V K16', nt_poin, vk16=v_sdappa_noms)
!
! - Datastructure for tangents at each node
!
    sdappa_tgno = sdappa(1:19)//'.TGNO'
    call wkvect(sdappa_tgno, 'V V R', 6*nb_cont_node, vr=v_sdappa_tgno)
!
! - Datastructure for tangents at each node by element
!
    sdappa_tgel = sdappa(1:19)//'.TGEL'
    call jecrec(sdappa_tgel, 'V V R', 'NU', 'CONTIG', 'VARIABLE', &
                nb_cont_elem)
    call jeecra(sdappa_tgel, 'LONT', 6*nt_elem_node)
    longt = 0
    do i_cont_elem = 1, nb_cont_elem
        elem_indx = i_cont_elem
        call cfnben(ds_contact%sdcont_defi, elem_indx, 'CONNEX', elem_nbnode)
        longc = 6*elem_nbnode
        call jeecra(jexnum(sdappa_tgel, i_cont_elem), 'LONMAX', ival=longc)
        call jecroc(jexnum(sdappa_tgel, i_cont_elem))
        longt = longt+longc
    end do
    ASSERT(longt .eq. 6*nt_elem_node)
!
! - Datastructure for check normals discontinuity
!
    sdappa_verk = sdappa(1:19)//'.VERK'
    sdappa_vera = sdappa(1:19)//'.VERA'
    call wkvect(sdappa_verk, 'V V K8', nb_node_mesh, vk8=v_sdappa_verk)
    call wkvect(sdappa_vera, 'V V R', nb_node_mesh, vr=v_sdappa_vera)
    call jeecra(sdappa_verk, 'LONUTI', 0)
!
    call jedema()
!
end subroutine
