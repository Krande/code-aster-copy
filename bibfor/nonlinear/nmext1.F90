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

subroutine nmext1(mesh, field, field_disc, field_type, field_s, &
                  nb_elem, nb_node, nb_poin, nb_spoi, nb_cmp, &
                  type_extr_elem, type_extr, type_extr_cmp, type_sele_cmp, &
                  list_node, list_elem, list_poin, list_spoi, list_cmp, &
                  work_node, work_poin, work_elem)
!
    implicit none
!
#include "asterfort/nmext2.h"
#include "asterfort/nmext3.h"
!
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1504
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: field
    character(len=4), intent(in) :: field_disc
    character(len=24), intent(in) :: field_type
    character(len=24), intent(in) :: field_s
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_poin
    integer(kind=8), intent(in) :: nb_spoi
    integer(kind=8), intent(in) :: nb_cmp
    character(len=8), intent(in) :: type_extr_elem
    character(len=8), intent(in) :: type_extr
    character(len=8), intent(in) :: type_extr_cmp
    character(len=8), intent(in) :: type_sele_cmp
    character(len=24), intent(in) :: list_node
    character(len=24), intent(in) :: list_elem
    character(len=24), intent(in) :: list_poin
    character(len=24), intent(in) :: list_spoi
    character(len=24), intent(in) :: list_cmp
    character(len=19), intent(in) :: work_node
    character(len=19), intent(in) :: work_poin
    character(len=19), intent(in) :: work_elem
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Extract values and store them in working vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  field            : name of field
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  field_type       : type of field (name in results datastructure)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  nb_elem          : number of elements
! In  nb_node          : number of nodes
! In  nb_poin          : number of points (Gauss)
! In  nb_spoi          : number of subpoints
! In  nb_cmp           : number of components
! In  type_extr_elem   : type of extraction by element
! In  type_extr        : type of extraction
! In  type_extr_cmp    : type of extraction for components
! In  type_sele_cmp    : type of selection for components NOM_CMP or NOM_VARI
! In  list_node        : name of object contains list of nodes
! In  list_elem        : name of object contains list of elements
! In  list_poin        : name of object contains list of points (Gauss)
! In  list_spoi        : name of object contains list of subpoints
! In  list_cmp         : name of object contains list of components
! In  work_node        : working vector to save node values
! In  work_poin        : working vector to save point (Gauss) values
! In  work_elem        : working vector to save element values
!
! --------------------------------------------------------------------------------------------------
!
!
! - For nodal values
!
    if (field_disc .eq. 'NOEU' .and. nb_node > 0) then
        call nmext2(mesh, field, nb_cmp, nb_node, type_extr, &
                    type_extr_cmp, list_node, list_cmp, work_node)
    end if
!
! - For point (Gauss) values
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        if (nb_elem > 0 .or. nb_poin > 0) then
            call nmext3(mesh, field, field_type, field_s, &
                        nb_cmp, nb_elem, nb_poin, nb_spoi, &
                        type_extr_elem, type_extr, type_extr_cmp, type_sele_cmp, &
                        list_elem, list_poin, list_spoi, list_cmp, &
                        work_poin, work_elem)
        end if
    end if
!
end subroutine
