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

subroutine dfc_read_disc(sdcont, zoneKeyword, mesh, model, model_ndim, &
                         nb_cont_zone, lLineRela, listRela)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfc_read_zone.h"
#include "asterfort/dfc_save_dime.h"
#include "asterfort/dfc_chck.h"
#include "asterfort/quadco.h"
#include "asterfort/elimcq.h"
#include "asterfort/tablco.h"
#include "asterfort/cacoco.h"
#include "asterfort/cacoeq.h"
#include "asterfort/capoco.h"
#include "asterfort/elimco.h"
#include "asterfort/sansco.h"
#include "asterfort/typeco.h"
!
    character(len=8), intent(in) :: sdcont, mesh, model
    character(len=16), intent(in) :: zoneKeyword
    integer(kind=8), intent(in) :: model_ndim, nb_cont_zone
    aster_logical, intent(out) :: lLineRela
    character(len=19), intent(out) :: listRela
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Discrete method - Read contact data
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  zoneKeyword      : factor keyword zone of contact
! In  mesh             : name of mesh
! In  model            : name of model
! In  model_ndim       : dimension of model
! In  nb_cont_zone     : number of zones of contact
! Out nb_cont_surf     : number of surfaces of contact
! Out nb_cont_elem     : number of elements of contact
! Out nb_cont_node     : number of nodes of contact
! Out lLineRela        : flag for linear relations
! Out listRela         : name of object for linear relations
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_cont_surf, nb_cont_elem, nb_cont_node
    aster_logical :: l_elim_coq3d, l_node_q8
!
! --------------------------------------------------------------------------------------------------
!
    l_node_q8 = .true.
    l_elim_coq3d = .true.
!
! - QUAD8 specific treatment activation (l_node_q8 = .true. if need linearization, see CACOEQ)
!
    call quadco(sdcont, l_node_q8)
!
! - Read zone: nodes and elements
!
    call dfc_read_zone(sdcont, zoneKeyword, mesh, model, nb_cont_zone, &
                       nb_cont_surf, nb_cont_elem, nb_cont_node)
!
! - Cleaning nodes and elements
!
    call elimco(sdcont, mesh, model, nb_cont_surf, &
                nb_cont_elem, nb_cont_node, l_elim_coq3d)
!
! - QUAD8 specific treatment: suppress middle nodes in contact lists
!
    if (l_node_q8) then
        call elimcq(sdcont, mesh, nb_cont_zone, nb_cont_surf, nb_cont_node)
    end if
!
! - Inverse connectivities
!
    call tablco(sdcont, mesh, nb_cont_surf, nb_cont_elem, nb_cont_node)
!
! - Save contact counters
!
    call dfc_save_dime(sdcont, mesh, model_ndim, nb_cont_zone, nb_cont_surf, &
                       nb_cont_elem, nb_cont_node)
!
! - Keyword SANS_GROUP_NO
!
    call sansco(sdcont, zoneKeyword, mesh)
!
! - Elements and nodes parameters
!
    call typeco(sdcont, mesh)
!
! - Some checks
!
    call dfc_chck(sdcont, mesh, model_ndim)
!
! - Gap for beams
!
    call capoco(sdcont, zoneKeyword)
!
! - Gap for shells
!
    call cacoco(sdcont, zoneKeyword, mesh)
!
! - Create QUAD8 linear relations
!
    if (l_node_q8) then
        call cacoeq(sdcont, mesh, lLineRela, listRela)
    end if
!
end subroutine
