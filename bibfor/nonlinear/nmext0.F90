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

subroutine nmext0(field_disc, nb_elem, nb_node, nb_poin, nb_spoi, &
                  nb_cmp, work_node, work_poin, work_elem, type_extr_elem, &
                  type_extr)
!
    implicit none
!
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=4), intent(in) :: field_disc
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_poin
    integer(kind=8), intent(in) :: nb_spoi
    integer(kind=8), intent(in) :: nb_cmp
    character(len=8), intent(in) :: type_extr_elem
    character(len=8), intent(in) :: type_extr
    character(len=19), intent(in) :: work_poin
    character(len=19), intent(in) :: work_node
    character(len=19), intent(in) :: work_elem
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Create workink vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  nb_node          : number of nodes
! In  nb_elem          : number of elements
! In  nb_poin          : number of points (Gauss)
! In  nb_spoi          : number of subpoints
! In  nb_cmp           : number of components
! In  work_node        : working vector to save node values
! In  work_elem        : working vector to save element values
! In  work_poin        : working vector to save point (Gauss) values
! In  type_extr        : type of extraction
! In  type_extr_elem   : type of extraction by element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: node_init_vale, elem_init_vale
    integer(kind=8) :: ino, i_elem, i_cmp, i_poin, i_spoi
    real(kind=8), pointer :: v_work_elem(:) => null()
    real(kind=8), pointer :: v_work_poin(:) => null()
    real(kind=8), pointer :: v_work_node(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    if (field_disc .eq. 'NOEU') then
        if (nb_node == 0) then
            go to 999
        end if
        call wkvect(work_node, 'V V R', nb_node*nb_cmp, vr=v_work_node)
    else if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        if (nb_poin == 0) then
            go to 999
        end if
        call wkvect(work_poin, 'V V R', nb_poin*nb_spoi*nb_cmp, vr=v_work_poin)
        if (nb_elem == 0) then
            go to 999
        end if
        call wkvect(work_elem, 'V V R', nb_elem*nb_poin*nb_spoi*nb_cmp, vr=v_work_elem)
    else
        ASSERT(.false.)
    end if
!
! - Values for initializations
!
    if (type_extr .eq. 'MAX') then
        node_init_vale = -r8maem()
    else if (type_extr .eq. 'MIN') then
        node_init_vale = +r8maem()
    else if (type_extr .eq. 'MAXI_ABS') then
        node_init_vale = 0.d0
    else if (type_extr .eq. 'MINI_ABS') then
        node_init_vale = +r8maem()
    else if (type_extr .eq. 'VALE') then
        node_init_vale = 0.d0
    else if (type_extr .eq. 'MOY') then
        node_init_vale = 0.d0
    else
        ASSERT(.false.)
    end if
!
    if (field_disc .eq. 'ELGA') then
        if (type_extr_elem .eq. 'MAX') then
            elem_init_vale = -r8maem()
        else if (type_extr_elem .eq. 'MIN') then
            elem_init_vale = +r8maem()
        else if (type_extr_elem .eq. 'VALE') then
            elem_init_vale = 0.d0
        else if (type_extr_elem .eq. 'MOY') then
            elem_init_vale = 0.d0
        else
            ASSERT(.false.)
        end if
    end if
!
    if (field_disc .eq. 'ELEM') then
        if (type_extr_elem .eq. 'VALE') then
            elem_init_vale = 0.d0
        else
            ASSERT(.false.)
        end if
    end if
!
! - Set for working vector (nodes)
!
    if (field_disc .eq. 'NOEU') then
        do ino = 1, nb_node
            do i_cmp = 1, nb_cmp
                v_work_node(nb_cmp*(ino-1)+i_cmp) = node_init_vale
            end do
        end do
    end if
!
! - Set for working vector (elements and points)
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        do i_elem = 1, nb_elem
            do i_poin = 1, nb_poin
                do i_spoi = 1, nb_spoi
                    do i_cmp = 1, nb_cmp
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem-1)+ &
                                    nb_poin*nb_spoi*(i_cmp-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = node_init_vale
                        v_work_poin(1+nb_poin*nb_spoi*(i_cmp-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = elem_init_vale
                    end do
                end do
            end do
        end do
    end if
!
999 continue
!
end subroutine
