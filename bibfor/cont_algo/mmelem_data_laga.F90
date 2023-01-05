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

subroutine mmelem_data_laga(l_axi_, &
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
    integer, parameter :: nb_cont_geom = 33
    integer, parameter :: nb_cont_solv = 33
!
! - Name of geometry type for slave element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_slav_name = (/ &
                                       'SEG2    ', 'SEG3    ', 'SEG2    ', 'SEG3    ', 'TRIA3   ', &
                                       'TRIA3   ', 'TRIA6   ', 'TRIA6   ', 'QUAD4   ', 'QUAD4   ', &
                                       'QUAD8   ', 'QUAD8   ', 'QUAD4   ', 'TRIA3   ', 'TRIA6   ', &
                                       'QUAD4   ', 'TRIA6   ', 'QUAD8   ', 'TRIA6   ', 'QUAD9   ', &
                                       'QUAD8   ', 'TRIA3   ', 'QUAD8   ', 'QUAD9   ', 'QUAD9   ', &
                                       'QUAD4   ', 'QUAD9   ', 'TRIA3   ', 'QUAD9   ', 'POI1    ', &
                                                            'POI1    ', 'POI1    ', 'POI1    '/)
!
! - Name of geometry type for master element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_mast_name = (/ &
                                       'SEG2    ', 'SEG3    ', 'SEG3    ', 'SEG2    ', 'TRIA3   ', &
                                       'TRIA6   ', 'TRIA3   ', 'TRIA6   ', 'QUAD4   ', 'QUAD8   ', &
                                       'QUAD4   ', 'QUAD8   ', 'TRIA3   ', 'QUAD4   ', 'QUAD4   ', &
                                       'TRIA6   ', 'QUAD8   ', 'TRIA6   ', 'QUAD9   ', 'TRIA6   ', &
                                       'TRIA3   ', 'QUAD8   ', 'QUAD9   ', 'QUAD8   ', 'QUAD4   ', &
                                       'QUAD9   ', 'TRIA3   ', 'QUAD9   ', 'QUAD9   ', 'LAG2    ', &
                                                            'NOLAG2  ', 'LAG3    ', 'NOLAG3  '/)
!
! - Name of geometry type for contact/friction element
!
    character(len=8), parameter, dimension(nb_cont_geom) :: lypg_cont_name = (/ &
                                       'SEG22   ', 'SEG33   ', 'SEG23   ', 'SEG32   ', 'TRIA33  ', &
                                       'TR3TR6  ', 'TR6TR3  ', 'TRIA66  ', 'QUAD44  ', 'QU4QU8  ', &
                                       'QU8QU4  ', 'QUAD88  ', 'QU4TR3  ', 'TR3QU4  ', 'TR6QU4  ', &
                                       'QU4TR6  ', 'TR6QU8  ', 'QU8TR6  ', 'TR6QU9  ', 'QU9TR6  ', &
                                       'QU8TR3  ', 'TR3QU8  ', 'QU8QU9  ', 'QU9QU8  ', 'QU9QU4  ', &
                                       'QU4QU9  ', 'QU9TR3  ', 'TR3QU9  ', 'QUAD99  ', 'POI1    ', &
                                                            'POI1    ', 'POI1    ', 'POI1    '/)
!
! - Number of nodes for contact/friction element
!
    integer, parameter, dimension(nb_cont_geom) :: nb_node = (/ &
                                                   4, 6, 5, 5, 6, &
                                                   9, 9, 12, 8, 12, &
                                                   12, 16, 7, 7, 10, &
                                                   10, 14, 14, 15, 15, &
                                                   11, 11, 17, 17, 13, &
                                                   13, 12, 12, 18, 1, &
                                                   1, 1, 1/)
!
! - Name of FE type of contact element
!
    character(len=8), parameter, dimension(nb_cont_solv) :: lypf_cont_name = (/ &
                                       'CMS2S2  ', 'CMS3S3  ', 'CMS2S3  ', 'CMS3S2  ', 'CMT3T3  ', &
                                       'CMT3T6  ', 'CMT6T3  ', 'CMT6T6  ', 'CMQ4Q4  ', 'CMQ4Q8  ', &
                                       'CMQ8Q4  ', 'CMQ8Q8  ', 'CMQ4T3  ', 'CMT3Q4  ', 'CMT6Q4  ', &
                                       'CMQ4T6  ', 'CMT6Q8  ', 'CMQ8T6  ', 'CMT6Q9  ', 'CMQ9T6  ', &
                                       'CMQ8T3  ', 'CMT3Q8  ', 'CMQ8Q9  ', 'CMQ9Q8  ', 'CMQ9Q4  ', &
                                       'CMQ4Q9  ', 'CMQ9T3  ', 'CMT3Q9  ', 'CMQ9Q9  ', 'CMP1L2  ', &
                                                            'CMP1N2  ', 'CMP1L3  ', 'CMP1N3  '/)
!
! - Name of FE type of friction element
!
    character(len=8), parameter, dimension(nb_cont_solv) :: lypf_frot_name = (/ &
                                       'FMS2S2  ', 'FMS3S3  ', 'FMS2S3  ', 'FMS3S2  ', 'FMT3T3  ', &
                                       'FMT3T6  ', 'FMT6T3  ', 'FMT6T6  ', 'FMQ4Q4  ', 'FMQ4Q8  ', &
                                       'FMQ8Q4  ', 'FMQ8Q8  ', 'FMQ4T3  ', 'FMT3Q4  ', 'FMT6Q4  ', &
                                       'FMQ4T6  ', 'FMT6Q8  ', 'FMQ8T6  ', 'FMT6Q9  ', 'FMQ9T6  ', &
                                       'FMQ8T3  ', 'FMT3Q8  ', 'FMQ8Q9  ', 'FMQ9Q8  ', 'FMQ9Q4  ', &
                                       'FMQ4Q9  ', 'FMQ9T3  ', 'FMT3Q9  ', 'FMQ9Q9  ', 'FMP1L2  ', &
                                                            'FMP1N2  ', 'FMP1L3  ', 'FMP1N3  '/)
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
