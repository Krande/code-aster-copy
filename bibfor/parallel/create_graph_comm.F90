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

subroutine create_graph_comm(object, nb_comm, comm, tag)
!
use sort_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/build_tree_comm.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: object
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
    character(len=16) :: type
    character(len=24) :: k24
    integer :: i, iret
    integer, pointer :: v_domdis(:) => null()
    integer, pointer :: v_comm(:) => null()
    integer, pointer :: v_tag(:) => null()
    integer, pointer :: v_joint(:) => null()
!
    call jemarq()
!
    call gettco(object, type)
!
    nb_comm = 0
!
! --- Result depends on type
    if(type == 'MAILLAGE_P') then
        k24 = object//'.DOMJOINTS'
        call jeexin(k24, iret)
        if( iret > 0 ) then
            call jelira(object//'.DOMJOINTS', 'LONMAX', nb_comm, k8bid)
            call jeveuo(object//'.DOMJOINTS', 'L', vi=v_domdis)
        end if
    else
        k24 = object//'.NUME.NBJO'
        call jeexin(k24, iret)
        if( iret > 0) then
            type = "NUME_DDL"
            call jeveuo(k24, 'L', vi=v_joint)
            nb_comm = v_joint(1)
            AS_ALLOCATE(vi=v_domdis, size=max(1, nb_comm))
            do i = 1, nb_comm
                v_domdis(i) = v_joint(i+1)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
! --- Allocation
    call wkvect(comm, 'V V I', max(1, nb_comm), vi=v_comm)
    call wkvect(tag, 'V V I', max(1, nb_comm), vi=v_tag)
!
! --- Create graph
    call build_tree_comm(v_domdis, nb_comm, v_comm, v_tag)
!
    if(type== 'NUME_DDL') then
        AS_DEALLOCATE(vi=v_domdis)
    end if
!
    call jedema()
!
end subroutine
