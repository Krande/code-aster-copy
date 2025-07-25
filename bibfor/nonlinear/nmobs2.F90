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

subroutine nmobs2(meshz, sd_obsv, tabl_name, time, title, &
                  field_disc, field_type, field_s, &
                  nb_elem, nb_node, nb_poin, nb_spoi, nb_cmp, &
                  type_extr_elem, type_extr, type_extr_cmp, type_sele_cmp, &
                  list_node, list_elem, list_poin, list_spoi, &
                  list_cmp, list_vari, &
                  field, work_node, work_elem, nb_obsf_effe)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmobsz.h"
#include "asterfort/sdmpic.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1504
!
    character(len=*), intent(in) :: meshz
    character(len=19), intent(in) :: sd_obsv
    character(len=19), intent(in) :: tabl_name
    real(kind=8), intent(in) :: time
    character(len=16), intent(in) :: title
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=24), intent(in) :: field_s
    character(len=4), intent(in) :: field_disc
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_poin
    integer(kind=8), intent(in) :: nb_spoi
    integer(kind=8), intent(in) :: nb_cmp
    character(len=24), intent(in) :: list_node
    character(len=24), intent(in) :: list_elem
    character(len=24), intent(in) :: list_poin
    character(len=24), intent(in) :: list_spoi
    character(len=24), intent(in) :: list_cmp
    character(len=24), intent(in) :: list_vari
    character(len=8), intent(in) :: type_extr
    character(len=8), intent(in) :: type_extr_elem
    character(len=8), intent(in) :: type_extr_cmp
    character(len=8), intent(in) :: type_sele_cmp
    character(len=19), intent(in) :: work_node
    character(len=19), intent(in) :: work_elem
    integer(kind=8), intent(inout) :: nb_obsf_effe
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear operators - Observation
!
! Compute values and save them in table
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  title            : title of observation
! In  sd_obsv          : datastructure for observation parameters
! In  tabl_name        : name of observation table
! In  time             : current time
! In  field            : name of field
! In  field_type       : type of field (name in results datastructure)
! In  field_disc       : localization of field (discretization: NOEU, ELGA or ELEM)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  list_node        : name of object contains list of nodes
! In  nb_node          : number of nodes
! In  list_elem        : name of object contains list of elements
! In  nb_elem          : number of elements
! In  list_poin        : name of object contains list of points (Gauss)
! In  nb_poin          : number of points (Gauss)
! In  list_spoi        : name of object contains list of subpoints
! In  nb_spoi          : number of subpoints
! In  list_cmp         : name of object contains list of components (NOM_CMP)
! In  list_vari        : name of object contains list of components (NOM_VARI)
! In  nb_cmp           : number of components
! In  type_extr        : type of extraction
! In  type_extr_elem   : type of extraction by element
! In  type_extr_cmp    : type of extraction for components
! In  type_sele_cmp    : type of selection for components NOM_CMP or NOM_VARI
! In  work_node        : working vector to save node values
! In  work_elem        : working vector to save element values
! IO  nb_obsf_effe     : number of _effective_observations
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) nb_para_maxi
    parameter(nb_para_maxi=20)
    character(len=16) :: v_cmp_name(nb_para_maxi)
!
    integer(kind=8) :: i_node, i_elem, i_poin, i_spoi, i_cmp
    integer(kind=8) :: iret
    real(kind=8) :: vale_r
    integer(kind=8) :: nb_node_r, nb_elem_r, nb_cmp_r, nb_poin_r, nb_spoi_r
    integer(kind=8) :: nb_poin_e, nb_spoi_e, nb_poin_elem, nb_spoi_elem
    integer(kind=8) :: poin_nume, spoi_nume, node_nume, elem_nume, nume_glob
    character(len=8) :: node_name, elem_name
    character(len=16) :: cmp_name
    aster_logical :: l_pmesh
    integer(kind=8), pointer :: cesd(:) => null()
    character(len=8), pointer :: v_list_cmp(:) => null()
    character(len=16), pointer :: v_list_vari(:) => null()
    integer(kind=8), pointer :: v_list_node(:) => null()
    integer(kind=8), pointer :: v_list_elem(:) => null()
    integer(kind=8), pointer :: v_list_poin(:) => null()
    integer(kind=8), pointer :: v_list_spoi(:) => null()
    integer(kind=8), pointer :: v_nonulg(:) => null()
    real(kind=8), pointer :: v_work_node(:) => null()
    real(kind=8), pointer :: v_work_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    l_pmesh = isParallelMesh(meshz)
!
! - Convert to reduced field
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        call jeexin(field_s, iret)
        if (iret .eq. 0) then
            call sdmpic('CHAM_ELEM', field)
            call celces(field, 'V', field_s)
        end if
        call jeveuo(field_s(1:19)//'.CESD', 'L', vi=cesd)
    end if
!
! - Number of nodes for loop
!
    if (field_disc .eq. 'NOEU') then
        if (type_extr .eq. 'VALE') then
            nb_node_r = nb_node
        elseif ((type_extr .eq. 'MIN') .or. &
                (type_extr .eq. 'MAX') .or. &
                (type_extr .eq. 'MAXI_ABS') .or. &
                (type_extr .eq. 'MINI_ABS') .or. &
                (type_extr .eq. 'MOY')) then
            nb_node_r = 1
        else
            ASSERT(.false.)
        end if
    end if
!
! - Number of elements for loop
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        if (type_extr .eq. 'VALE') then
            nb_elem_r = nb_elem
        elseif ((type_extr .eq. 'MIN') .or. &
                (type_extr .eq. 'MAX') .or. &
                (type_extr .eq. 'MAXI_ABS') .or. &
                (type_extr .eq. 'MINI_ABS') .or. &
                (type_extr .eq. 'MOY')) then
            nb_elem_r = 1
        else
            ASSERT(.false.)
        end if
    end if
!
! - Number for components for loop
!
    if (type_extr_cmp .eq. ' ') then
        nb_cmp_r = nb_cmp
    else
        nb_cmp_r = 1
    end if
!
! - Get name of components
!
    call jeveuo(list_cmp, 'L', vk8=v_list_cmp)
    if (type_sele_cmp .eq. 'NOM_VARI') then
        call jeveuo(list_vari, 'L', vk16=v_list_vari)
    end if
    ASSERT(nb_cmp .le. nb_para_maxi)
    do i_cmp = 1, nb_cmp
        if (type_sele_cmp .eq. 'NOM_CMP') then
            v_cmp_name(i_cmp) = v_list_cmp(i_cmp)
        elseif (type_sele_cmp .eq. 'NOM_VARI') then
            v_cmp_name(i_cmp) = v_list_vari(i_cmp)
        else
            ASSERT(.false.)
        end if
    end do
!
! - For node discretization
!
    if (field_disc .eq. 'NOEU' .and. nb_node > 0) then
        call jeveuo(work_node, 'L', vr=v_work_node)
        call jeveuo(list_node, 'L', vi=v_list_node)
        if (l_pmesh) then
            call jeveuo(meshz(1:8)//'.NUNOLG', 'L', vi=v_nonulg)
        end if
!
        do i_node = 1, nb_node_r
!
! --------- Current node
!
            node_nume = v_list_node(i_node)
            node_name = int_to_char8(node_nume)
            if (l_pmesh) then
                nume_glob = v_nonulg(node_nume)
            end if
!
! --------- Write values
!
            do i_cmp = 1, nb_cmp_r
                vale_r = v_work_node(i_cmp+nb_cmp*(i_node-1))
                cmp_name = v_cmp_name(i_cmp)
                if (l_pmesh) then

                    call nmobsz(sd_obsv, tabl_name, title, field_type, field_disc, &
                                type_extr, type_extr_cmp, type_extr_elem, type_sele_cmp, cmp_name, &
                                time, vale_r, &
                                node_namez=node_name, glob_numez=nume_glob)
                else
                    call nmobsz(sd_obsv, tabl_name, title, field_type, field_disc, &
                                type_extr, type_extr_cmp, type_extr_elem, type_sele_cmp, cmp_name, &
                                time, vale_r, &
                                node_namez=node_name)
                end if
                nb_obsf_effe = nb_obsf_effe+1
            end do
        end do
    end if
!
! - For element discretization
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        if (nb_elem > 0) then
            call jeveuo(work_elem, 'L', vr=v_work_elem)
            call jeveuo(list_elem, 'L', vi=v_list_elem)
        end if
!
        if (nb_poin > 0) then
            call jeveuo(list_poin, 'L', vi=v_list_poin)
            call jeveuo(list_spoi, 'L', vi=v_list_spoi)
        end if
!
        do i_elem = 1, nb_elem_r

            if (type_extr .eq. 'VALE') then
!
! ------------- Current element
!
                elem_nume = v_list_elem(i_elem)
                elem_name = int_to_char8(elem_nume)
!
! ------------- Real number of point/subpoint for current element
!
                nb_poin_elem = cesd(1+5+4*(elem_nume-1))
                nb_spoi_elem = cesd(1+5+4*(elem_nume-1)+1)
!
! ------------- Check
!
                nb_poin_e = nb_poin
                nb_spoi_e = nb_spoi
                if (nb_poin_e .gt. nb_poin_elem) nb_poin_e = nb_poin_elem
                if (nb_spoi_e .gt. nb_spoi_elem) nb_spoi_e = nb_spoi_elem
!
! ------------- Number for points/subpoints for loop
!
                if (type_extr_elem .eq. 'VALE') then
                    nb_poin_r = nb_poin_e
                    nb_spoi_r = nb_spoi_e
                else
                    nb_poin_r = 1
                    nb_spoi_r = 1
                end if
            else
                nb_poin_r = 1
                nb_spoi_r = 1
            end if
!
            do i_poin = 1, nb_poin_r
                do i_spoi = 1, nb_spoi_r
!
! ----------------- Current number of point/subpoint
!
                    if (type_extr_elem .eq. 'VALE') then
                        poin_nume = v_list_poin(i_poin)
                        spoi_nume = v_list_spoi(i_spoi)
                    else
                        poin_nume = i_poin
                        spoi_nume = i_spoi
                    end if
!
! ----------------- Write values
!
                    do i_cmp = 1, nb_cmp_r
                        vale_r = v_work_elem(nb_cmp*nb_poin*nb_spoi*(i_elem-1)+ &
                                             nb_poin*nb_spoi*(i_cmp-1)+ &
                                             nb_spoi*(i_poin-1)+ &
                                             (i_spoi-1)+1)
                        cmp_name = v_cmp_name(i_cmp)
                        call nmobsz(sd_obsv, tabl_name, title, field_type, field_disc, &
                                    type_extr, type_extr_cmp, type_extr_elem, &
                                    type_sele_cmp, cmp_name, &
                                    time, vale_r, &
                                    elem_namez=elem_name, &
                                    poin_numez=poin_nume, spoi_numez=spoi_nume)
                        nb_obsf_effe = nb_obsf_effe+1
                    end do
                end do
            end do
        end do
    end if
!
end subroutine
