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

subroutine cfleqc(mesh, sdcont_defi, nb_cont_zone, nb_cont_node, nb_cont_surf, &
                  v_poin_node, v_indi_node, nb_node_elim)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/cfnbsf.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mminfl.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sdcont_defi
    integer(kind=8), intent(in) :: nb_cont_zone
    integer(kind=8), intent(in) :: nb_cont_surf
    integer(kind=8), intent(in) :: nb_cont_node
    integer(kind=8), pointer :: v_poin_node(:)
    integer(kind=8), pointer :: v_indi_node(:)
    integer(kind=8), intent(out) :: nb_node_elim
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Suppress quadratic middle nodes of QUAD8 - Create list of (middle) nodes to suppress
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  nb_cont_zone     : number of zones of contact
! In  nb_cont_surf     : number of surfaces of contact
! In  nb_cont_node     : number of nodes of contact (after detection of middle nodes)
! Out v_indi_node      : pointer to indicator of middle nodes
! Out v_poin_node      : pointer to pointer of contact surface
! Out nb_node_elim     : number of nodes to suppress
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: elem_nume, type_nume
    integer(kind=8) :: i_surf_curr, i_zone, i_elem, i_node, i_node_quad, i_surf
    integer(kind=8) :: nb_surf, nb_elem, nb_node_quad, nb_node
    character(len=8) :: type_name
    aster_logical :: l_veri
    integer(kind=8) :: jdecma, jdecno, jdecqu
    integer(kind=8) :: node_nume_1, node_nume_2
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=24) :: sdcont_noeuco
    integer(kind=8), pointer :: v_sdcont_noeuco(:) => null()
    character(len=24) :: sdcont_pzoneco
    integer(kind=8), pointer :: v_sdcont_pzoneco(:) => null()
    character(len=24) :: sdcont_pnoeuqu
    integer(kind=8), pointer :: v_sdcont_pnoeuqu(:) => null()
    character(len=24) :: sdcont_noeuqu
    integer(kind=8), pointer :: v_sdcont_noeuqu(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_node_elim = 0
    jdecqu = 0
!
! - Create vectors to use in cfmeno subroutine
!
    AS_ALLOCATE(vi=v_poin_node, size=nb_cont_surf+1)
    AS_ALLOCATE(vi=v_indi_node, size=nb_cont_node)
!
! - Datastructure for contact definition
!
    sdcont_pzoneco = sdcont_defi(1:16)//'.PZONECO'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    sdcont_noeuco = sdcont_defi(1:16)//'.NOEUCO'
    sdcont_pnoeuqu = sdcont_defi(1:16)//'.PNOEUQU'
    sdcont_noeuqu = sdcont_defi(1:16)//'.NOEUQU'
    call jeveuo(sdcont_pzoneco, 'L', vi=v_sdcont_pzoneco)
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
    call jeveuo(sdcont_noeuco, 'L', vi=v_sdcont_noeuco)
    call jeveuo(sdcont_pnoeuqu, 'L', vi=v_sdcont_pnoeuqu)
    call jeveuo(sdcont_noeuqu, 'L', vi=v_sdcont_noeuqu)
!
! - Access to mesh
!
    call jeveuo(mesh(1:8)//'.TYPMAIL', 'L', vi=v_mesh_typmail)
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- No computation
!
        l_veri = mminfl(sdcont_defi, 'VERIF', i_zone)
        if (l_veri) then
            goto 21
        end if
!
! ----- Number of contact surfaces
!
        nb_surf = v_sdcont_pzoneco(i_zone+1)-v_sdcont_pzoneco(i_zone)
        ASSERT(nb_surf .eq. 2)
!
! ----- Loop on surfaces
!
        do i_surf = 1, nb_surf
!
! --------- Parameters of current surface
!
            i_surf_curr = nb_surf*(i_zone-1)+i_surf
            call cfnbsf(sdcont_defi, i_surf_curr, 'MAIL', nb_elem, jdecma)
!
! --------- Change pointer
!
            v_poin_node(i_surf_curr+1) = v_poin_node(i_surf_curr)
!
! --------- Loop on elements
!
            do i_elem = 1, nb_elem
!
! ------------- Current element
!
                elem_nume = v_sdcont_mailco(jdecma+i_elem)
!
! ------------- Type of element
!
                type_nume = v_mesh_typmail(elem_nume)
                call jenuno(jexnum('&CATA.TM.NOMTM', type_nume), type_name)
!
! ------------- Suppress middle nodes of QUAD8
!
                if (type_name(1:5) .eq. 'QUAD8') then
                    nb_node_quad = (v_sdcont_pnoeuqu(i_zone+1)-v_sdcont_pnoeuqu(i_zone))/3
                    jdecqu = v_sdcont_pnoeuqu(i_zone)
                    call cfnbsf(sdcont_defi, i_surf_curr, 'NOEU', nb_node, jdecno)
                    do i_node_quad = 1, nb_node_quad
                        node_nume_1 = v_sdcont_noeuqu(jdecqu+3*(i_node_quad-1)+1)
                        do i_node = 1, nb_node
                            node_nume_2 = v_sdcont_noeuco(jdecno+i_node)
                            if (node_nume_1 .eq. node_nume_2) then
                                if (v_indi_node(jdecno+i_node) .eq. 0) then
                                    v_indi_node(jdecno+i_node) = 1
                                    v_poin_node(i_surf_curr+1) = v_poin_node(i_surf_curr+1)+1
                                    nb_node_elim = nb_node_elim+1
                                end if
                            end if
                        end do
                    end do
                end if
            end do
        end do
21      continue
    end do
!
    ASSERT((2*nb_cont_zone) .eq. nb_cont_surf)
!
end subroutine
