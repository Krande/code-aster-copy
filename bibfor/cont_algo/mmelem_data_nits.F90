! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine mmelem_data_nits(l_axi_, &
                            typg_slav_name_, typg_mast_name_, &
                            nb_cont_type_, nb_node_elem_, &
                            typg_cont_nume_, &
                            typf_cont_nume_, &
                            typf_frot_nume_, &
                            get_elem_indx_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
!
    aster_logical, intent(in) :: l_axi_
    character(len=8), intent(in) :: typg_slav_name_
    character(len=8), intent(in) :: typg_mast_name_
    integer, intent(out) :: nb_cont_type_
    integer, intent(out) :: nb_node_elem_
    integer, intent(out) :: typg_cont_nume_
    integer, intent(out) :: typf_cont_nume_
    integer, intent(out) :: typf_frot_nume_
    integer, intent(out) :: get_elem_indx_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Mortar method - Define contact/friction elements
!
! --------------------------------------------------------------------------------------------------
!
! In  l_axi            : .true. for axi-symetric model
! In  typg_slav_name   : name of geometric type of slave element
! In  typg_mast_name   : name of geometric type of master element
! Out nb_cont_type     : total number of contact elements defined
! Out nb_node_elem     : number of nodes of contact/friction element
! Out typg_cont_nume   : index of geometric type of contact/friction element
! Out typf_cont_nume   : index of FE type of contact element
! Out typf_frot_nume   : index of FE type of friction element
! Out get_elem_indx    : index to select contact/friction element (get)
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nb_cont_geom = 10
    integer, parameter :: nb_cont_solv = 10
!
! - Name of geometry type for slave element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_slav_name = (/ &
                                                   'TRIA3   ', 'TRIA3   ', 'TRIA6   ', 'TRIA6   ', &
                                                   'QUAD4   ', 'QUAD4   ', 'QUAD8   ', 'QUAD8   ', &
                                                            'QUAD9   ', 'QUAD9   '/)
!
! - Name of geometry type for master element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_mast_name = (/ &
                                                   'SEG2    ', 'SEG3    ', 'SEG2    ', 'SEG3    ', &
                                                   'SEG2    ', 'SEG3    ', 'SEG2    ', 'SEG3    ', &
                                                            'SEG2    ', 'SEG3    '/)
!
! - Name of geometry type for contact/friction element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_cont_name = (/ &
                                                   'TR3SE2  ', 'TR3SE3  ', 'TR6SE2  ', 'TR6SE3  ', &
                                                   'QU4SE2  ', 'QU4SE3  ', 'QU8SE2  ', 'QU8SE3  ', &
                                                            'QU9SE2  ', 'QU9SE3  '/)
!
! - Number of nodes for contact/friction element
!
    integer, parameter, dimension(nb_cont_geom) :: nb_node = (/ &
                                                   5, 6, 8, 9, &
                                                   6, 7, 10, 11, &
                                                   11, 12/)
!
! - Name of FE type of contact element
!
    character(len=8), parameter, dimension(nb_cont_solv) :: lypf_cont_name = (/ &
                                                   'CNT3S2  ', 'CNT3S3  ', 'CNT6S2  ', 'CNT6S3  ', &
                                                   'CNQ4S2  ', 'CNQ4S3  ', 'CNQ8S2  ', 'CNQ8S3  ', &
                                                            'CNQ9S2  ', 'CMQ9S3  '/)
!
! - Name of FE type of friction element
!
    character(len=8), parameter, dimension(nb_cont_solv) :: lypf_frot_name = (/ &
                                                   'FNT3S2  ', 'FNT3S3  ', 'FNT6S2  ', 'FNT6S3  ', &
                                                   'FNQ4S2  ', 'FNQ4S3  ', 'FNQ8S2  ', 'FNQ8S3  ', &
                                                            'FNQ9S2  ', 'FMQ9S3  '/)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i_cont_geom, elem_indx
    character(len=16) :: typg_cont_name, valk(2), typf_cont_name, typf_frot_name
!
! --------------------------------------------------------------------------------------------------
!
    elem_indx = 0
!
! - Total number of contact elements defined
!
    nb_cont_type_ = nb_cont_solv
!
! - Index to select contact/friction element
!
    elem_indx = 0
    do i_cont_geom = 1, nb_cont_geom
        if (typg_slav_name_ .eq. lypg_slav_name(i_cont_geom)) then
            if (typg_mast_name_ .eq. lypg_mast_name(i_cont_geom)) then
                elem_indx = i_cont_geom
            end if
        end if
    end do
    if (elem_indx .eq. 0) then
        valk(1) = typg_slav_name_
        valk(2) = typg_mast_name_
        call utmess('F', 'CONTACT_96', nk=2, valk=valk)
    end if
!
! - Number of nodes of contact/friction element
!
    ASSERT(elem_indx .ne. 0)
    nb_node_elem_ = nb_node(elem_indx)
!
! - Index of geometric type of contact/friction element
!
    ASSERT(elem_indx .ne. 0)
    typg_cont_name = lypg_cont_name(elem_indx)
    call jenonu(jexnom('&CATA.TM.NOMTM', typg_cont_name), typg_cont_nume_)
!
! - Name of FE type of contact element
!
    ASSERT(elem_indx .ne. 0)
    typf_cont_name = lypf_cont_name(elem_indx)
    if (l_axi_) then
        typf_cont_name(7:7) = 'A'
    end if
    call jenonu(jexnom('&CATA.TE.NOMTE', typf_cont_name), typf_cont_nume_)
!
! - Name of FE type of friction element
!
    ASSERT(elem_indx .ne. 0)
    typf_frot_name = lypf_frot_name(elem_indx)
    if (l_axi_) then
        typf_frot_name(7:7) = 'A'
    end if
    call jenonu(jexnom('&CATA.TE.NOMTE', typf_frot_name), typf_frot_nume_)
!
! - Get index for contact/friction element
!
    ASSERT(elem_indx .ne. 0)
    get_elem_indx_ = elem_indx
!
end subroutine
