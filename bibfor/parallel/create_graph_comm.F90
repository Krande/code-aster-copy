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

subroutine create_graph_comm(object, type, nb_comm, comm, tag)
!
    use sort_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/build_tree_comm.h"
#include "asterfort/dismoi.h"
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
    integer(kind=8), intent(out) :: nb_comm
    character(len=*), intent(in) :: comm, tag
!
!---------------------------------------------------------------------------------------------------
!
! Le but est de construire le graphe de comm optimisé pour les comm point à point
!
!---------------------------------------------------------------------------------------------------
!
    character(len=8) :: k8bid
    character(len=19) :: k19
    character(len=24) :: k24, pgid, gcom
    integer(kind=8) :: iret
    integer(kind=8), pointer :: v_domdis(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
    integer(kind=8), pointer :: v_gcom(:) => null()
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    mpi_int :: mpicou
!
    call jemarq()
!
    nb_comm = 0
!
! --- Result depends on type
    if (type == 'MAILLAGE_P') then
        k24 = object(1:8)//'.JOIN      .DOMJ'
        pgid = object(1:8)//'.JOIN      .PGID'
        gcom = object(1:8)//'.JOIN      .GCOM'
    elseif (type == "NUME_EQUA") then
        call dismoi("JOINTS", object, "NUME_EQUA", repk=k19, arret="F")
        k24 = k19//".DOMJ"
        pgid = k19//".PGID"
        gcom = k19//".GCOM"
    elseif (type == "LIGREL") then
        call dismoi("JOINTS", object, "LIGREL", repk=k19, arret="F")
        k24 = k19//".DOMJ"
        pgid = k19//".PGID"
        gcom = k19//".GCOM"
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeexin(k24, iret)
    if (iret > 0) then
        call jelira(k24, 'LONUTI', nb_comm, k8bid)
        call jeveuo(k24, 'L', vi=v_domdis)
        call jeveuo(pgid, 'L', vi4=v_pgid)
        call jeveuo(gcom, 'L', vi=v_gcom)
        mpicou = v_gcom(1)
    else
        mpicou = -1
    end if
!
! --- Allocation
    call wkvect(comm, 'V V I', max(1, nb_comm), vi=v_comm)
    call wkvect(tag, 'V V I', max(1, nb_comm), vi=v_tag)
!
! --- Create graph
    call build_tree_comm(v_domdis, nb_comm, v_pgid, mpicou, v_comm, v_tag)
!
    call jedema()
!
end subroutine
