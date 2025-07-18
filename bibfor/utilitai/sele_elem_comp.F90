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

subroutine sele_elem_comp(modelz, compor, defo_comp, list_elem_comp)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/as_allocate.h"
#include "asterfort/exisdg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/Behaviour_type.h"
!
!
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: compor
    character(len=16), intent(in) :: defo_comp
    integer(kind=8), pointer :: list_elem_comp(:)
!
! --------------------------------------------------------------------------------------------------
!
! Select elements by comportment
!
! --------------------------------------------------------------------------------------------------
!
! In  modelz          : name of model
! In  compor          : name of comportment CARTE
! In  defo_comp       : type of deformation to find
! In  list_elem_comp  : elements preselected on complete mesh
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mesh
    character(len=16) :: grandeur_name, defo_comp_f
    character(len=19) :: ligrmo
    character(len=24) :: name_liel
    integer(kind=8) :: nume_elem
    integer(kind=8) :: idx_gd, idx_cmp, nb_gd_max, iret, nb_cmp_max, nb_elem_mesh
    integer(kind=8) :: nb_elem_grel, nb_grel, nb_ec
    integer(kind=8) :: i_elem_grel, i_grel
    integer(kind=8), pointer :: comp_ptma(:) => null()
    integer(kind=8), pointer :: comp_desc(:) => null()
    character(len=16), pointer :: comp_vale(:) => null()
    integer(kind=8), pointer :: list_elem_grel(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    ligrmo = modelz(1:8)//'.MODELE'
    call dismoi('NOM_MAILLA', modelz, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nb_elem_mesh)
!
! - Allocate list of elements
!
    AS_ALLOCATE(vi=list_elem_comp, size=nb_elem_mesh)
!
! - Preparation of comportment datas
!
    grandeur_name = 'COMPOR'
    call dismoi('NB_EC', grandeur_name, 'GRANDEUR', repi=nb_ec)
    ASSERT(nb_ec .le. 1)
    call etenca(compor, ligrmo, iret)
    ASSERT(iret .eq. 0)
    call jeveuo(compor(1:19)//'.DESC', 'L', vi=comp_desc)
    nb_gd_max = comp_desc(2)
    call jelira(jexnom('&CATA.GD.NOMCMP', grandeur_name), 'LONMAX', nb_cmp_max)
    call jeveuo(compor(1:19)//'.VALE', 'L', vk16=comp_vale)
    call jeveuo(compor(1:19)//'.PTMA', 'L', vi=comp_ptma)
!
! - Select elements in list
!
    nb_grel = nbgrel(ligrmo)
    name_liel = ligrmo//'.LIEL'
    do i_grel = 1, nb_grel
        nb_elem_grel = nbelem(ligrmo, i_grel)
        call jeveuo(jexnum(name_liel, i_grel), 'L', vi=list_elem_grel)
        do i_elem_grel = 1, nb_elem_grel
            nume_elem = list_elem_grel(i_elem_grel)
            idx_gd = comp_ptma(nume_elem)
            if (idx_gd .ne. 0) then
                idx_cmp = comp_desc(3+2*nb_gd_max+idx_gd)
                ASSERT(exisdg([idx_cmp], 1))
                if (exisdg([idx_cmp], 1)) then
                    defo_comp_f = comp_vale((idx_gd-1)*nb_cmp_max+DEFO)
                    if (defo_comp .eq. defo_comp_f) then
                        list_elem_comp(nume_elem) = 1
                    end if
                end if
            end if
        end do
    end do
!
end subroutine
