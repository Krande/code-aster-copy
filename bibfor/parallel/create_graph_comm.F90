! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine create_graph_comm(object, type, nb_comm, comm, tag)
!
use sort_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/build_tree_comm.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: object
    character(len=*), intent(in) :: type
    integer, intent(inout) :: nb_comm
    character(len=*), intent(in) :: comm, tag
!
!---------------------------------------------------------------------------------------------------
!
! Le but est de construire le graphe de comm optimisé pour les comm point à point
!
!---------------------------------------------------------------------------------------------------
!
    character(len=8) :: k8bid
    character(len=24) :: k24
    integer :: iret
    integer, pointer :: v_domdis(:) => null()
    integer, pointer :: v_comm(:) => null()
    integer, pointer :: v_tag(:) => null()
!
    call jemarq()
!
    nb_comm = 0
!
! --- Result depends on type
    if(type == 'MAILLAGE_P') then
        k24 = object//'.DOMJOINTS'
    elseif(type == "NUME_DDL") then
        k24 = object//'.NUME.DOMJ'
    elseif(type == "LIGREL") then
        k24 = object//'.DOMJ'
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeexin(k24, iret)
    if( iret > 0 ) then
        call jelira(k24, 'LONUTI', nb_comm, k8bid)
        call jeveuo(k24, 'L', vi=v_domdis)
    end if
!
! --- Allocation
    call wkvect(comm, 'V V I', max(1, nb_comm), vi=v_comm)
    call wkvect(tag, 'V V I', max(1, nb_comm), vi=v_tag)
!
! --- Create graph
    call build_tree_comm(v_domdis, nb_comm, v_comm, v_tag)
!
    call jedema()
!
end subroutine
