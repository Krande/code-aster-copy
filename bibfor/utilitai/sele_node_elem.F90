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

subroutine sele_node_elem(modelz, nb_elem_type, list_elem_type, list_node, nb_node_found, &
                          pre_select_elem)
!
    implicit none
!
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/typele.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/dismoi.h"
!
!
    character(len=*), intent(in) :: modelz
    integer(kind=8), intent(in) :: nb_elem_type
    character(len=16), pointer :: list_elem_type(:)
    integer(kind=8), pointer :: list_node(:)
    integer(kind=8), intent(out) :: nb_node_found
    integer(kind=8), pointer, optional :: pre_select_elem(:)
!
! --------------------------------------------------------------------------------------------------
!
! Select nodes by element type
!
! --------------------------------------------------------------------------------------------------
!
! In  modelz          : name of model
! In  nb_elem_type    : number of element to detect
! In  list_elem_type  : list of element to detect
! Out list_node       : list of nodes where element type has been detected
! Out nb_node_found   : number of nodes where element type has been detected
! In  pre_select_elem : elements preselected on complete mesh
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mesh
    character(len=16) :: type_elem
    character(len=19) :: ligrmo
    character(len=24) :: name_liel
    integer(kind=8) :: nb_elem_grel, nb_grel, nb_node_elem, nb_node_mesh
    integer(kind=8) :: i_grel, i_node, i_elem_grel, i_node_elem, i_elem
    integer(kind=8) :: idx_type_elem, nume_node, nume_elem
    logical :: l_find, l_sele
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: list_elem_grel(:) => null()
    integer(kind=8), pointer :: list_node_all(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_node_found = 0
    ligrmo = modelz(1:8)//'.MODELE'
    call dismoi('NOM_MAILLA', modelz, 'MODELE', repk=mesh)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node_mesh)
!
! - Allocate list of nodes
!
    AS_ALLOCATE(vi=list_node_all, size=nb_node_mesh)
    AS_ALLOCATE(vi=list_node, size=nb_node_mesh)
!
! - Create list of nodes with special elements
!
    nb_grel = nbgrel(ligrmo)
    name_liel = ligrmo//'.LIEL'
    do i_grel = 1, nb_grel
        idx_type_elem = typele(ligrmo, i_grel)
        call jenuno(jexnum('&CATA.TE.NOMTE', idx_type_elem), type_elem)
        l_find = .false.
        do i_elem = 1, nb_elem_type
            if (type_elem .eq. list_elem_type(i_elem)) then
                l_find = .true.
                goto 10
            end if
        end do
10      continue
        if (l_find) then
            nb_elem_grel = nbelem(ligrmo, i_grel)
            call jeveuo(jexnum(name_liel, i_grel), 'L', vi=list_elem_grel)
            do i_elem_grel = 1, nb_elem_grel
                nume_elem = list_elem_grel(i_elem_grel)
                if (present(pre_select_elem)) then
                    l_sele = .false.
                    if ((pre_select_elem(nume_elem)) .ne. 0) then
                        l_sele = .true.
                    end if
                else
                    l_sele = .true.
                end if
                if (l_sele) then
                    call jeveuo(jexnum(mesh//'.CONNEX', nume_elem), 'L', vi=connex)
                    call jelira(jexnum(mesh//'.CONNEX', nume_elem), 'LONMAX', nb_node_elem)
                    do i_node_elem = 1, nb_node_elem
                        nume_node = connex(i_node_elem)
                        list_node_all(nume_node) = 1
                    end do
                end if
            end do
        end if
    end do
!
    i_node = 0
    do nume_node = 1, nb_node_mesh
        if (list_node_all(nume_node) .eq. 1) then
            i_node = i_node+1
            list_node(i_node) = nume_node
        end if
    end do
    nb_node_found = i_node
!
    AS_DEALLOCATE(vi=list_node_all)
!
end subroutine
