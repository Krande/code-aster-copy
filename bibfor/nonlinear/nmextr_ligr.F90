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
subroutine nmextr_ligr(meshz, modelz, sdextrz, nb_keyw_fact, nb_field_comp)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/jeveuo.h"
!
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: meshz
    character(len=*), intent(in) :: sdextrz
    integer(kind=8), intent(in) :: nb_keyw_fact
    integer(kind=8), intent(in) :: nb_field_comp
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Field extraction datastructure
!
! Create LIGREL for fields not a default in nonlinear operator
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  sdextr           : name of datastructure for extraction
! In  nb_keyw_fact     : number of factor keyword to read extraction parameters
! In  nb_field_comp    : number of fields to compute (not a default in nonlinear operator)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_keyw_fact, i_field_comp, i_elem, i_elem_mesh, i_elem_ligr
    integer(kind=8) :: nb_elem_mesh, nb_elem_ligr, nb_elem, nume_elem
    character(len=2) :: chaine
    character(len=24) :: ligrel, list_elem
    character(len=24) :: field_disc, field_type, field_comp, field
    character(len=14) :: sdextr
    character(len=24) :: extr_comp, extr_info, extr_field
    character(len=24), pointer :: v_extr_comp(:) => null()
    character(len=24), pointer :: v_extr_field(:) => null()
    integer(kind=8), pointer :: v_extr_info(:) => null()
    integer(kind=8), pointer :: list_elem_mesh(:) => null()
    integer(kind=8), pointer :: list_elem_ligr(:) => null()
    integer(kind=8), pointer :: v_list_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdextr = sdextrz
    call dismoi('NB_MA_MAILLA', meshz, 'MAILLAGE', repi=nb_elem_mesh)
!
! - Access to datastructure
!
    extr_info = sdextr(1:14)//'     .INFO'
    call jeveuo(extr_info, 'L', vi=v_extr_info)
    extr_comp = sdextr(1:14)//'     .COMP'
    call jeveuo(extr_comp, 'L', vk24=v_extr_comp)
    extr_field = sdextr(1:14)//'     .CHAM'
    call jeveuo(extr_field, 'L', vk24=v_extr_field)
!
! - Create list of elements
!
    AS_ALLOCATE(vi=list_elem_mesh, size=nb_elem_mesh)
!
    do i_field_comp = 1, nb_field_comp
!
        list_elem_mesh(1:nb_elem_mesh) = 0
!
! ----- Get info
!
        field_comp = v_extr_comp(4*(i_field_comp-1)+1)
        field_disc = v_extr_comp(4*(i_field_comp-1)+2)
        field_type = v_extr_comp(4*(i_field_comp-1)+3)
        ligrel = v_extr_comp(4*(i_field_comp-1)+4)
        ASSERT(field_type .eq. 'EPSI_ELGA')
        ASSERT(field_disc .eq. 'ELGA')
!
! ----- Loop on all observations
!
        nb_elem_ligr = 0
        do i_keyw_fact = 1, nb_keyw_fact
            field = v_extr_field(4*(i_keyw_fact-1)+4)
            if (field .eq. field_comp) then
                write (chaine, '(I2)') i_keyw_fact
                list_elem = sdextr(1:14)//chaine(1:2)//'   .MAIL'
                call jeveuo(list_elem, 'L', vi=v_list_elem)
                nb_elem = v_extr_info(7+7*(i_keyw_fact-1)+3)
                do i_elem = 1, nb_elem
                    nume_elem = v_list_elem(i_elem)
                    if (list_elem_mesh(nume_elem) .eq. 0) then
                        nb_elem_ligr = nb_elem_ligr+1
                        list_elem_mesh(nume_elem) = 1
                    end if
                end do
            end if
        end do
!
! ----- Create list of elements
!
        AS_ALLOCATE(vi=list_elem_ligr, size=nb_elem_ligr)
        i_elem_ligr = 0
        do i_elem_mesh = 1, nb_elem_mesh
            nume_elem = i_elem_mesh
            if (list_elem_mesh(nume_elem) .eq. 1) then
                i_elem_ligr = i_elem_ligr+1
                list_elem_ligr(i_elem_ligr) = nume_elem
            end if
        end do
        ASSERT(i_elem_ligr .eq. nb_elem_ligr)
!
! ----- Create LIGREL
!
        call exlim1(list_elem_ligr, nb_elem_ligr, modelz, 'V', ligrel)
        AS_DEALLOCATE(vi=list_elem_ligr)
    end do
!
    AS_DEALLOCATE(vi=list_elem_mesh)
!
end subroutine
