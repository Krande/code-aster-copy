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
!
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine addPhantomNodesFromCells(mesh, indic_nodes)
!
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i4.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/create_graph_comm.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
!
    character(len=8), intent(in) :: mesh
    integer(kind=4), intent(inout) :: indic_nodes(*)
!
! ---------------------------------------------------------------------------------------------
!
! Le but est de partager l'information sur les noeuds fantômes
! pour avoir la liste exhaustive des noeuds appartements à des GROUP_MA non-locaux
!
! ---------------------------------------------------------------------------------------------
!
    character(len=8) :: k8bid
    character(len=4) :: chnbjo
    character(len=24) :: nojoie, nojoir
    character(len=19) :: tag_name, comm_name
    integer :: rang, nbproc, nb_comm, domdis, iret
    integer :: i_comm, nbnoee, nbnoer, i_no, numno
    mpi_int :: mrank, msize, mpicou, tag4, numpr4, n4e, n4r
    integer, pointer :: v_comm(:) => null()
    integer, pointer :: v_tag(:) => null()
    integer, pointer :: v_joine(:) => null()
    integer, pointer :: v_joinr(:) => null()
    integer(kind=4), pointer :: v_send(:) => null()
    integer(kind=4), pointer :: v_recv(:) => null()
!
    call jemarq()
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    ASSERT(nbproc <= MT_DOMMAX)
    DEBUG_MPI('addPhantomNodesFromCells', rang, nbproc)
!
! --- Lecture des joints
    call jeexin(mesh//'.DOMJOINTS', iret)
    if(iret > 0) then
        comm_name = '&&ADDNODES.COMM'
        tag_name = '&&ADDNODES.TAG'
        call create_graph_comm(mesh, "MAILLAGE_P", nb_comm, comm_name, tag_name)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
!
        do i_comm = 1, nb_comm
            domdis = v_comm(i_comm)
! --- Get JOINT
            call codlet(domdis, 'G', chnbjo)
            nojoie = mesh//'.E'//chnbjo
            nojoir = mesh//'.R'//chnbjo
            call jeveuo(nojoie, 'L', vi=v_joine)
            call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
            call jeveuo(nojoir, 'L', vi=v_joinr)
            call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
            nbnoee = nbnoee/2
            nbnoer = nbnoer/2
!
            call wkvect('&&APNFCS.NOSEND', 'V V S', nbnoee, vi4=v_send)
            call wkvect('&&APNFCS.NORECV', 'V V S', nbnoer, vi4=v_recv)
!
! --- Get nodes to send
!
            do i_no = 1, nbnoee
                numno = v_joine(2*(i_no-1)+1)
                v_send(i_no) = indic_nodes(numno)
            end do
!
! --- Send and Recive
            n4e = to_mpi_int(nbnoee)
            n4r = to_mpi_int(nbnoer)
            tag4 = to_mpi_int(v_tag(i_comm))
            numpr4 = to_mpi_int(domdis)
            call asmpi_sendrecv_i4(v_send, n4e, numpr4, tag4, &
                                   v_recv, n4r, numpr4, tag4, mpicou)
!
! --- Write nodes to receive
!
            do i_no = 1, nbnoer
                numno = v_joinr(2*(i_no-1)+1)
                indic_nodes(numno) = max(indic_nodes(numno), v_recv(i_no))
            end do
!
! --- Cleaning
            call jedetr('&&APNFCS.NOSEND')
            call jedetr('&&APNFCS.NORECV')
        end do
!
! --- Cleaning
        call jedetr(comm_name)
        call jedetr(tag_name)
!
    endif
!
    call jedema()
!
end subroutine
