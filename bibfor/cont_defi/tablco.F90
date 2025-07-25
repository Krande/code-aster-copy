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

subroutine tablco(sdcont, mesh, nb_cont_surf, nb_cont_elem, nb_cont_node)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/juveca.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: sdcont
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nb_cont_surf
    integer(kind=8), intent(in) :: nb_cont_elem
    integer(kind=8), intent(in) :: nb_cont_node
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Inverse connectivity
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  mesh             : name of mesh
! In  nb_cont_surf     : number of surfaces of contact
! In  nb_cont_node     : number of nodes of contact
! In  nb_cont_elem     : number of elements of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_elem, i_node, i_cont_surf, i
    integer(kind=8) :: node_nume, elem_nume, numalo, no
    integer(kind=8) :: nb_node, nb_elem, nb_node_elem
    integer(kind=8) :: jdecno, jdecma
    integer(kind=8) :: pre_length, inc, long
    integer(kind=8) :: nt_elem_node, nt_node_elem
    integer(kind=8), pointer :: v_work_node(:) => null()
    integer(kind=8), pointer :: v_work_elem(:) => null()
    integer(kind=8), pointer :: v_connex(:) => null()
    integer(kind=8), pointer :: v_connex_longcum(:) => null()
    character(len=19) :: connex_inv
    integer(kind=8), pointer :: v_coninv(:) => null()
    integer(kind=8), pointer :: v_coninv_longcum(:) => null()
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=24) :: sdcont_noeuco
    integer(kind=8), pointer :: v_sdcont_noeuco(:) => null()
    character(len=24) :: sdcont_pzoneco
    integer(kind=8), pointer :: v_sdcont_pzoneco(:) => null()
    character(len=24) :: sdcont_psunoco
    integer(kind=8), pointer :: v_sdcont_psunoco(:) => null()
    character(len=24) :: sdcont_psumaco
    integer(kind=8), pointer :: v_sdcont_psumaco(:) => null()
    character(len=24) :: sdcont_manoco
    integer(kind=8), pointer :: v_sdcont_manoco(:) => null()
    character(len=24) :: sdcont_pmanoco
    integer(kind=8), pointer :: v_sdcont_pmanoco(:) => null()
    character(len=24) :: sdcont_nomaco
    integer(kind=8), pointer :: v_sdcont_nomaco(:) => null()
    character(len=24) :: sdcont_pnomaco
    integer(kind=8), pointer :: v_sdcont_pnomaco(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
! Access to number of nodes of an element
#define elem_nbnode(i_elem) v_connex_longcum(i_elem+1)-v_connex_longcum(i_elem)
!
! Access to connectivity of element
#define numglm(i_elem,i_node) v_connex(v_connex_longcum(i_elem)+i_node-1)
!
! --------------------------------------------------------------------------------------------------
!
    pre_length = 20*max(nb_cont_node, nb_cont_elem)
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    sdcont_noeuco = sdcont_defi(1:16)//'.NOEUCO'
    sdcont_pzoneco = sdcont_defi(1:16)//'.PZONECO'
    sdcont_psumaco = sdcont_defi(1:16)//'.PSUMACO'
    sdcont_psunoco = sdcont_defi(1:16)//'.PSUNOCO'
    sdcont_manoco = sdcont_defi(1:16)//'.MANOCO'
    sdcont_pmanoco = sdcont_defi(1:16)//'.PMANOCO'
    sdcont_nomaco = sdcont_defi(1:16)//'.NOMACO'
    sdcont_pnomaco = sdcont_defi(1:16)//'.PNOMACO'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
    call jeveuo(sdcont_noeuco, 'L', vi=v_sdcont_noeuco)
    call jeveuo(sdcont_pzoneco, 'L', vi=v_sdcont_pzoneco)
    call jeveuo(sdcont_psumaco, 'E', vi=v_sdcont_psumaco)
    call jeveuo(sdcont_psunoco, 'E', vi=v_sdcont_psunoco)
!
! - Access to mesh
!
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_longcum)
!
! - Construct inverse connectivity
!
    connex_inv = '&&TABLCO.CONINV'
    call cncinv(mesh, v_sdcont_mailco, nb_cont_elem, 'V', connex_inv)
    call jeveuo(connex_inv, 'L', vi=v_coninv)
    call jeveuo(jexatr(connex_inv, 'LONCUM'), 'L', vi=v_coninv_longcum)
!
! - Working vectors
!
    AS_ALLOCATE(vi=v_work_node, size=nb_cont_node)
    AS_ALLOCATE(vi=v_work_elem, size=nb_cont_elem)
!
! - Construct MANOCO objects (list of elements connected to each contact node)
!
    call wkvect(sdcont_manoco, 'G V I', pre_length, vi=v_sdcont_manoco)
    call wkvect(sdcont_pmanoco, 'G V I', nb_cont_node+1, vi=v_sdcont_pmanoco)
    v_sdcont_pmanoco(1) = 0
!
! - Loop on surfaces
!
    inc = 0
    long = pre_length
    do i_cont_surf = 1, nb_cont_surf
!
! ----- Pointers to nodes/elements of current surface
!
        jdecno = v_sdcont_psunoco(i_cont_surf)
        jdecma = v_sdcont_psumaco(i_cont_surf)
!
! ----- Numbers of nodes/elements of current surface
!
        nb_node = v_sdcont_psunoco(i_cont_surf+1)-v_sdcont_psunoco(i_cont_surf)
        nb_elem = v_sdcont_psumaco(i_cont_surf+1)-v_sdcont_psumaco(i_cont_surf)
!
! ----- Loop on nodes
!
        do i_node = 1, nb_node
            v_work_node(jdecno+i_node) = 0
!
! --------- Current node
!
            node_nume = v_sdcont_noeuco(jdecno+i_node)
!
! --------- Loop on elements connected to current node
!
            nb_node_elem = v_coninv_longcum(node_nume+1)-v_coninv_longcum(node_nume)
            do i_elem = 1, nb_node_elem
                numalo = v_coninv(v_coninv_longcum(node_nume)+i_elem-1)
                if (numalo .le. v_sdcont_psumaco(i_cont_surf+1) .and. &
                    numalo .gt. v_sdcont_psumaco(i_cont_surf)) then
                    inc = inc+1
                    if (inc .gt. long) then
                        long = 2*long
                        call juveca(sdcont_manoco, long)
                        call jeveuo(sdcont_manoco, 'E', vi=v_sdcont_manoco)
                    end if
                    v_sdcont_manoco(inc) = numalo
                    v_work_node(jdecno+i_node) = v_work_node(jdecno+i_node)+1
                end if
            end do
!
! --------- Update pointer
!
            v_sdcont_pmanoco(jdecno+i_node+1) = v_sdcont_pmanoco(jdecno+i_node)+ &
                                                v_work_node(jdecno+i_node)
        end do
    end do
!
! - Update length
!
    nt_elem_node = v_sdcont_pmanoco(nb_cont_node+1)
    ASSERT(nt_elem_node .le. long)
    call jeecra(sdcont_manoco, 'LONUTI', nt_elem_node)
!
! - Construct NOMACO objects (list of nodes connected to each contact element)
!
    call wkvect(sdcont_pnomaco, 'G V I', nb_cont_elem+1, vi=v_sdcont_pnomaco)
    v_sdcont_pnomaco(1) = 0
!
! - Compute length of NOMACO - Loop on surfaces
!
    inc = 0
    long = pre_length
    do i_cont_surf = 1, nb_cont_surf
!
! ----- Pointer to nodes of current surface
!
        jdecno = v_sdcont_psunoco(i_cont_surf)
!
! ----- Number of nodes of current surface
!
        nb_node = v_sdcont_psunoco(i_cont_surf+1)-v_sdcont_psunoco(i_cont_surf)
!
! ----- Loop on nodes
!
        do i_node = 1, nb_node
!
! --------- Current node
!
            node_nume = v_sdcont_noeuco(jdecno+i_node)
            nb_node_elem = v_coninv_longcum(node_nume+1)-v_coninv_longcum(node_nume)
!
! --------- Loop on elements connected to current node
!
            do i_elem = 1, nb_node_elem
                numalo = v_coninv(v_coninv_longcum(node_nume)+i_elem-1)
                if (numalo .le. v_sdcont_psumaco(i_cont_surf+1) .and. &
                    numalo .gt. v_sdcont_psumaco(i_cont_surf)) then
                    inc = inc+1
                    if (inc .gt. long) then
                        long = 2*long
                    end if
                end if
            end do
        end do
    end do
!
! - Upadte NOMACO - Loop on surfaces
!
    call wkvect(sdcont_nomaco, 'G V I', long, vi=v_sdcont_nomaco)
    inc = 0
    do i_cont_surf = 1, nb_cont_surf
!
! ----- Pointers to nodes/elements of current surface
!
        jdecno = v_sdcont_psunoco(i_cont_surf)
        jdecma = v_sdcont_psumaco(i_cont_surf)
!
! ----- Numbers of nodes/elements of current surface
!
        nb_node = v_sdcont_psunoco(i_cont_surf+1)-v_sdcont_psunoco(i_cont_surf)
        nb_elem = v_sdcont_psumaco(i_cont_surf+1)-v_sdcont_psumaco(i_cont_surf)
!
! ----- Loop on nodes
!
        do i_node = 1, nb_node
!
! --------- Current node
!
            node_nume = v_sdcont_noeuco(jdecno+i_node)
            nb_node_elem = v_coninv_longcum(node_nume+1)-v_coninv_longcum(node_nume)
!
! --------- Loop on elements connected to current node
!
            do i_elem = 1, nb_node_elem
                numalo = v_coninv(v_coninv_longcum(node_nume)+i_elem-1)
                if (numalo .le. v_sdcont_psumaco(i_cont_surf+1) .and. &
                    numalo .gt. v_sdcont_psumaco(i_cont_surf)) then
                    v_work_elem(numalo) = v_work_elem(numalo)+1
                end if
            end do
        end do
!
! ----- Loop on elements
!
        do i_elem = 1, nb_elem
!
! --------- Update pointer
!
            v_sdcont_pnomaco(jdecma+i_elem+1) = v_sdcont_pnomaco(jdecma+i_elem)+ &
                                                v_work_elem(jdecma+i_elem)
!
! --------- Current element
!
            elem_nume = v_sdcont_mailco(jdecma+i_elem)
!
! --------- Save nodes
!
            do i = 1, elem_nbnode(elem_nume)
                no = numglm(elem_nume, i)
                do i_node = 1, nb_node
                    node_nume = v_sdcont_noeuco(jdecno+i_node)
                    if (no .eq. node_nume) then
                        inc = inc+1
                        v_sdcont_nomaco(inc) = jdecno+i_node
                        goto 130
                    end if
                end do
130             continue
            end do
        end do
    end do
!
! - Update length
!
    nt_node_elem = v_sdcont_pnomaco(nb_cont_elem+1)
    ASSERT(nt_node_elem .le. long)
    call jeecra(sdcont_nomaco, 'LONUTI', nt_node_elem)
!
! - Clean
!
    AS_DEALLOCATE(vi=v_work_node)
    AS_DEALLOCATE(vi=v_work_elem)
    call jedetr(connex_inv)
!
end subroutine
