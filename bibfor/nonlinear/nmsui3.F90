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

subroutine nmsui3(ds_print, field_disc, nb_elem, nb_node, nb_poin, &
                  nb_spoi, nb_cmp, type_extr, type_extr_cmp, type_extr_elem, &
                  list_elem, work_node, work_elem, field, field_s, &
                  i_dof_monitor)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmsuiy.h"
#include "asterfort/sdmpic.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_poin
    integer(kind=8), intent(in) :: nb_spoi
    integer(kind=8), intent(in) :: nb_cmp
    character(len=24), intent(in) :: list_elem
    character(len=19), intent(in) :: field
    character(len=4), intent(in) :: field_disc
    character(len=24), intent(in) :: field_s
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=8), intent(in) :: type_extr
    character(len=8), intent(in) :: type_extr_elem
    character(len=8), intent(in) :: type_extr_cmp
    character(len=19), intent(in) :: work_node
    character(len=19), intent(in) :: work_elem
    integer(kind=8), intent(inout) :: i_dof_monitor
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear operators - DOF monitor
!
! Print monitored values in table
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_print         : datastructure for printing parameters
! In  nb_node          : number of nodes
! In  nb_elem          : number of elements
! In  nb_poin          : number of points (Gauss)
! In  nb_spoi          : number of subpoints
! In  nb_cmp           : number of components
! In  list_elem        : name of object contains list of elements
! In  field            : name of field
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  type_extr        : type of extraction
! In  type_extr_elem   : type of extraction by element
! In  type_extr_cmp    : type of extraction for components
! In  work_node        : working vector to save node values
! In  work_elem        : working vector to save element values
! IO  i_dof_monitor    : index of current monitoring
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_node, i_elem, i_poin, i_spoi, i_cmp
    integer(kind=8) :: nb_poin_elem, nb_spoi_elem, elem_nume, iret
    real(kind=8) :: vale_r
    integer(kind=8) :: nb_cmp_r, nb_poin_r, nb_spoi_r, nb_node_r, nb_elem_r
    integer(kind=8) :: nb_poin_e, nb_spoi_e
    integer(kind=8), pointer :: cesd(:) => null()
    integer(kind=8), pointer :: v_list_elem(:) => null()
    real(kind=8), pointer :: v_work_node(:) => null()
    real(kind=8), pointer :: v_work_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
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
! - For node discretization
!
    if (field_disc .eq. 'NOEU') then
        call jeveuo(work_node, 'L', vr=v_work_node)
        do i_node = 1, nb_node_r
            do i_cmp = 1, nb_cmp_r
                vale_r = v_work_node(i_cmp+nb_cmp*(i_node-1))
                call nmsuiy(ds_print, vale_r, i_dof_monitor)
            end do
        end do
    end if
!
! - For element discretization
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        call jeveuo(work_elem, 'L', vr=v_work_elem)
        call jeveuo(list_elem, 'L', vi=v_list_elem)
!
        do i_elem = 1, nb_elem_r
!
! --------- Current element
!
            elem_nume = v_list_elem(i_elem)
!
! --------- Real number of point/subpoint for current element
!
            nb_poin_elem = cesd(1+5+4*(elem_nume-1))
            nb_spoi_elem = cesd(1+5+4*(elem_nume-1)+1)
!
! --------- Check
!
            nb_poin_e = nb_poin
            nb_spoi_e = nb_spoi
            if (nb_poin_e .gt. nb_poin_elem) nb_poin_e = nb_poin_elem
            if (nb_spoi_e .gt. nb_spoi_elem) nb_spoi_e = nb_spoi_elem
!
! --------- Number for points/subpoints for loop
!
            if (type_extr_elem .eq. 'VALE') then
                nb_poin_r = nb_poin_e
                nb_spoi_r = nb_spoi_e
            else
                nb_poin_r = 1
                nb_spoi_r = 1
            end if
!
            do i_poin = 1, nb_poin_r
                do i_spoi = 1, nb_spoi_r
                    do i_cmp = 1, nb_cmp_r
                        vale_r = v_work_elem(nb_cmp*nb_poin*nb_spoi*(i_elem-1)+ &
                                             nb_poin*nb_spoi*(i_cmp-1)+ &
                                             nb_spoi*(i_poin-1)+ &
                                             (i_spoi-1)+1)
                        call nmsuiy(ds_print, vale_r, i_dof_monitor)
                    end do
                end do
            end do
        end do
    end if
!
end subroutine
