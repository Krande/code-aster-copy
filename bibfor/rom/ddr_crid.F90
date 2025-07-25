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
subroutine ddr_crid(ds_para, nb_node_rid, v_node_rid)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/addGroupElem.h"
#include "asterfort/addGroupNode.h"
!
    type(ROM_DS_ParaDDR), intent(in) :: ds_para
    integer(kind=8), intent(in)           :: nb_node_rid
    integer(kind=8), intent(in)           :: v_node_rid(nb_node_rid)
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_DOMAINE_REDUIT - Main process
!
! Construction of RID in mesh
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_para          : datastructure for parameters of EIM computation
! In  nb_node_rid      : number of nodes in RID
! In  v_node_rid       : list of nodes in RID
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_elem, nb_node, node_nbelem, elem_nbnode
    integer(kind=8) :: nunolo
    integer(kind=8) :: i_layer
    integer(kind=8) :: i_elem, i_node, i_elem_node, i_node_elem, i_rid_maxi
    integer(kind=8) :: nb_rid_elem, nb_int_node, nb_group_add, nb_sub_node, nb_layer_sub, i_couche
    integer(kind=8) :: indx, node_nume, elem_nume, nb_rid_node
    integer(kind=8) :: nb_layer_rid
    character(len=8) :: mesh
    character(len=24):: grelem_rid, grnode_int, grnode_sub
    integer(kind=8), pointer :: v_coninv(:) => null()
    integer(kind=8), pointer :: v_coninv_longcum(:) => null()
    integer(kind=8), pointer :: v_connex(:) => null()
    integer(kind=8), pointer :: v_connex_longcum(:) => null()
    aster_logical :: test, l_corr_ef
    aster_logical, pointer :: v_liel_rid(:) => null()
    aster_logical, pointer :: v_lino_rid(:) => null()
    aster_logical, pointer :: v_lino_add(:) => null()
    aster_logical, pointer :: v_list_in(:) => null()
    aster_logical, pointer :: v_list_sb(:) => null()
    aster_logical, pointer :: v_loca_sb(:) => null()
    integer(kind=8), pointer :: v_group(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM4_21')
    end if
!
! - Get parameters
!
    mesh = ds_para%mesh
    nb_layer_rid = ds_para%nb_layer_rid
    grelem_rid = ds_para%grelem_rid
    grnode_int = ds_para%grnode_int
    l_corr_ef = ds_para%l_corr_ef
    grnode_sub = ds_para%grnode_sub
    nb_layer_sub = ds_para%nb_layer_sub
!
! - Initializations
!
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', nb_elem)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', nb_node)
!
! - Access to mesh
!
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_longcum)
!
! - Access to inverse connectivity
!
    call cncinv(mesh, [0], 0, 'V', '&&MESHMA')
    call jeveuo('&&MESHMA', 'L', vi=v_coninv)
    call jeveuo(jexatr('&&MESHMA', 'LONCUM'), 'L', vi=v_coninv_longcum)
!
! - Create working objects
!
    AS_ALLOCATE(vl=v_liel_rid, size=nb_elem)
    AS_ALLOCATE(vl=v_lino_rid, size=nb_node)
    AS_ALLOCATE(vl=v_lino_add, size=nb_node)
    AS_ALLOCATE(vl=v_list_in, size=nb_node)
    AS_ALLOCATE(vl=v_list_sb, size=nb_node)
    AS_ALLOCATE(vl=v_loca_sb, size=nb_node)
!
! - Add GROUP_MA in mesh
!
    call addGroupElem(mesh, 1)
!
! - Add GROUP_NO in mesh
!
    if (ds_para%l_corr_ef) then
        nb_group_add = 2
    else
        nb_group_add = 1
    end if
    call addGroupNode(mesh, nb_group_add)
!
! - Set list of which nodes from mesh are in RID
!
    do i_node = 1, nb_node_rid
        node_nume = v_node_rid(i_node)
        ASSERT(node_nume .gt. 0)
        v_lino_rid(node_nume) = ASTER_TRUE
        v_lino_add(node_nume) = ASTER_TRUE
    end do
!
! - Add new nodes in RID from NB_LAYER
!
    do i_layer = 1, nb_layer_rid+1
        do i_node = 1, nb_node
            if (v_lino_rid(i_node)) then
                node_nbelem = v_coninv_longcum(i_node+1)-v_coninv_longcum(i_node)
                do i_node_elem = 1, node_nbelem
                    elem_nume = v_coninv(v_coninv_longcum(i_node)+i_node_elem-1)
                    elem_nbnode = v_connex_longcum(elem_nume+1)-v_connex_longcum(elem_nume)
                    do i_elem_node = 1, elem_nbnode
                        nunolo = v_connex(v_connex_longcum(elem_nume)+i_elem_node-1)
                        v_lino_add(nunolo) = ASTER_TRUE
                    end do
                end do
            end if
        end do
        do i_node = 1, nb_node
            v_lino_rid(i_node) = v_lino_add(i_node)
        end do
    end do
!
! - If all nodes of same element is in RID => this element is in RID
!
    do i_elem = 1, nb_elem
        elem_nbnode = v_connex_longcum(i_elem+1)-v_connex_longcum(i_elem)
        test = ASTER_FALSE
        do i_elem_node = 1, elem_nbnode
            nunolo = v_connex(v_connex_longcum(i_elem)+i_elem_node-1)
            if (v_lino_rid(nunolo)) then
                test = ASTER_TRUE
            else
                test = ASTER_FALSE
                exit
            end if
        end do
        if (test) then
            v_liel_rid(i_elem) = ASTER_TRUE
        end if
    end do
!
! - In case of existing limit domain: if element is in DOMAINE_MAXI => is in RID
!
    if (ds_para%l_rid_maxi) then
        do i_elem = 1, nb_elem
            if (v_liel_rid(i_elem)) then
                do i_rid_maxi = 1, ds_para%nb_rid_maxi
                    if (i_elem .eq. ds_para%v_rid_maxi(i_rid_maxi)) then
                        v_liel_rid(i_elem) = ASTER_TRUE
                        exit
                    else
                        v_liel_rid(i_elem) = ASTER_FALSE
                    end if
                end do
            end if
        end do
    end if
!
! - In case of existing limit domain: recreate the list of nodes in RID
!
    if (ds_para%l_rid_maxi) then
        AS_DEALLOCATE(vl=v_lino_rid)
        AS_ALLOCATE(vl=v_lino_rid, size=nb_node)
        do i_elem = 1, nb_elem
            if (v_liel_rid(i_elem)) then
                elem_nbnode = v_connex_longcum(i_elem+1)-v_connex_longcum(i_elem)
                do i_elem_node = 1, elem_nbnode
                    nunolo = v_connex(v_connex_longcum(i_elem)+i_elem_node-1)
                    v_lino_rid(nunolo) = ASTER_TRUE
                end do
            end if
        end do
    end if
!
! - Number of nodes/elements in RID
!
    nb_rid_elem = 0
    do i_elem = 1, nb_elem
        if (v_liel_rid(i_elem)) then
            nb_rid_elem = nb_rid_elem+1
        end if
    end do
    nb_rid_node = 0
    do i_node = 1, nb_node
        if (v_lino_rid(i_node)) then
            nb_rid_node = nb_rid_node+1
        end if
    end do
    if (niv .ge. 2) then
        call utmess('I', 'ROM4_22', si=nb_rid_elem)
        call utmess('I', 'ROM4_29', si=nb_rid_node)
    end if
!
! - Create group for elements of RID
!
    call jecroc(jexnom(mesh//'.GROUPEMA', grelem_rid))
    call jeecra(jexnom(mesh//'.GROUPEMA', grelem_rid), 'LONMAX', max(1, nb_rid_elem))
    call jeecra(jexnom(mesh//'.GROUPEMA', grelem_rid), 'LONUTI', nb_rid_elem)
    call jeveuo(jexnom(mesh//'.GROUPEMA', grelem_rid), 'E', vi=v_group)
    indx = 0
    do i_elem = 1, nb_elem
        if (v_liel_rid(i_elem)) then
            indx = indx+1
            v_group(indx) = i_elem
        end if
    end do
!
! - Create list of nodes of interface (in case of existing limit domain or not?)
!
    if (ds_para%l_rid_maxi) then
        do i_rid_maxi = 1, ds_para%nb_rid_maxi
            i_elem = ds_para%v_rid_maxi(i_rid_maxi)
            if (.not. v_liel_rid(i_elem)) then
                elem_nbnode = v_connex_longcum(i_elem+1)-v_connex_longcum(i_elem)
                do i_elem_node = 1, elem_nbnode
                    nunolo = v_connex(v_connex_longcum(i_elem)+i_elem_node-1)
                    if (v_lino_rid(nunolo)) then
                        v_list_in(nunolo) = .true._1
                        v_list_sb(nunolo) = .true._1
                        v_loca_sb(nunolo) = .true._1
                    end if
                end do
            end if
        end do
    else
        do i_elem = 1, nb_elem
            if (.not. v_liel_rid(i_elem)) then
                elem_nbnode = v_connex_longcum(i_elem+1)-v_connex_longcum(i_elem)
                do i_elem_node = 1, elem_nbnode
                    nunolo = v_connex(v_connex_longcum(i_elem)+i_elem_node-1)
                    if (v_lino_rid(nunolo)) then
                        v_list_in(nunolo) = .true._1
                        v_list_sb(nunolo) = .true._1
                        v_loca_sb(nunolo) = .true._1
                    end if
                end do
            end if
        end do
    end if
!
! - Number of nodes of interface
!
    nb_int_node = 0
    do i_node = 1, nb_node
        if (v_list_in(i_node)) then
            nb_int_node = nb_int_node+1
        end if
    end do
    if (niv .ge. 2) then
        call utmess('I', 'ROM4_23', si=nb_int_node)
    end if
!
! - Create group for nodes of interface
!
    call jecroc(jexnom(mesh//'.GROUPENO', grnode_int))
    call jeecra(jexnom(mesh//'.GROUPENO', grnode_int), 'LONMAX', max(1, nb_int_node))
    call jeecra(jexnom(mesh//'.GROUPENO', grnode_int), 'LONUTI', nb_int_node)
    call jeveuo(jexnom(mesh//'.GROUPENO', grnode_int), 'E', vi=v_group)
    if (nb_int_node .eq. 0) then
        call utmess('A', 'ROM4_27')
    end if
    indx = 0
    do i_node = 1, nb_node
        if (v_list_in(i_node)) then
            indx = indx+1
            v_group(indx) = i_node
        end if
    end do
!
! - Create list of nodes for EF correctors
!
    if (l_corr_ef) then
        do i_couche = 1, nb_layer_sub+1
            do i_node = 1, nb_node
                if (v_list_sb(i_node)) then
                    node_nbelem = v_coninv_longcum(i_node+1)-v_coninv_longcum(i_node)
                    do i_node_elem = 1, node_nbelem
                        elem_nume = v_coninv(v_coninv_longcum(i_node)+i_node_elem-1)
                        if (v_liel_rid(elem_nume)) then
                            elem_nbnode = v_connex_longcum(elem_nume+1)-v_connex_longcum(elem_nume)
                            do i_elem_node = 1, elem_nbnode
                                nunolo = v_connex(v_connex_longcum(elem_nume)+i_elem_node-1)
                                v_loca_sb(nunolo) = .true._1
                            end do
                        end if
                    end do
                end if
            end do
            do i_node = 1, nb_node
                v_list_sb(i_node) = v_loca_sb(i_node)
            end do
        end do
    end if
!
! - Number of nodes for EF correcteurs
!
    nb_sub_node = 0
    do i_node = 1, nb_node
        if (v_list_sb(i_node)) then
            nb_sub_node = nb_sub_node+1
        end if
    end do
    if (ds_para%l_corr_ef .and. nb_sub_node .eq. 0) then
        call utmess('A', 'ROM4_28')
    end if
    if (niv .ge. 2 .and. l_corr_ef) then
        call utmess('I', 'ROM4_26', si=nb_sub_node)
    end if
!
! - Create group for nodes of rigid zone in EF correctors
!
    if (l_corr_ef) then
        call jecroc(jexnom(mesh//'.GROUPENO', grnode_sub))
        call jeecra(jexnom(mesh//'.GROUPENO', grnode_sub), 'LONMAX', max(1, nb_sub_node))
        call jeecra(jexnom(mesh//'.GROUPENO', grnode_sub), 'LONUTI', nb_sub_node)
        call jeveuo(jexnom(mesh//'.GROUPENO', grnode_sub), 'E', vi=v_group)
        indx = 0
        do i_node = 1, nb_node
            if (v_list_sb(i_node)) then
                indx = indx+1
                v_group(indx) = i_node
            end if
        end do
    end if
!
! - Clean
!
    AS_DEALLOCATE(vl=v_liel_rid)
    AS_DEALLOCATE(vl=v_lino_rid)
    AS_DEALLOCATE(vl=v_lino_add)
    AS_DEALLOCATE(vl=v_list_in)
    AS_DEALLOCATE(vl=v_list_sb)
    AS_DEALLOCATE(vl=v_loca_sb)
!
end subroutine
