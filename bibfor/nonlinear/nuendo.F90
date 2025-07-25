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
subroutine nuendo(modelz, nume_ddl, sdnuen)
!
    implicit none
!
#include "asterfort/dismoi.h"
#include "asterfort/sele_node_elem.h"
#include "asterfort/select_dof.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: nume_ddl
    character(len=24), intent(in) :: sdnuen
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Initializations
!
! Get position of damaged dof
!
! --------------------------------------------------------------------------------------------------
!
! In  modelz   : name of model
! In  nume_ddl : name of numbering (NUME_DDL)
! In  sdnuen   : name of datastructure to save position of damaged dof
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_elem_type, nb_node_found, nbEqua, nb_cmp
    integer(kind=8), pointer :: listNode(:) => null()
    character(len=8), pointer :: listCmp(:) => null()
    integer(kind=8), pointer :: listEqua(:) => null()
    character(len=16), pointer :: list_elem_type(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
!
! - Create list of elements type
!
    nb_elem_type = 7
    AS_ALLOCATE(vk16=list_elem_type, size=nb_elem_type)
    list_elem_type(1) = 'MNDPTR6'
    list_elem_type(2) = 'MNDPQS8'
    list_elem_type(3) = 'MNAXTR6'
    list_elem_type(4) = 'MNAXQS8'
    list_elem_type(5) = 'MNVG_HEXA20'
    list_elem_type(6) = 'MNVG_TETRA10'
    list_elem_type(7) = 'MNVG_PENTA15'
!
! - Create list of components
!
    nb_cmp = 1
    AS_ALLOCATE(vk8=listCmp, size=nb_cmp)
    listCmp(1) = 'DAMG'
!
! - Select nodes by element type
!
    call sele_node_elem(modelz, nb_elem_type, list_elem_type, listNode, nb_node_found)
!
! - Create list of equations
!
    call dismoi('NB_EQUA', nume_ddl, 'NUME_DDL', repi=nbEqua)
    if (nb_node_found .gt. 0) then
        call wkvect(sdnuen, 'V V I', nbEqua, vi=listEqua)
    else
        goto 999
    end if
!
! - Find components in list of equations
!
    call select_dof(listEqua, &
                    numeDofZ_=nume_ddl, &
                    nbNodeToSelect_=nb_node_found, listNodeToSelect_=listNode, &
                    nbCmpToSelect_=nb_cmp, listCmpToSelect_=listCmp)
!
999 continue
!
    AS_DEALLOCATE(vi=listNode)
    AS_DEALLOCATE(vk8=listCmp)
    AS_DEALLOCATE(vk16=list_elem_type)
!
end subroutine
