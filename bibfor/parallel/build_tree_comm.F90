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

subroutine build_tree_comm(domdist, nbdom, comm, tag)
!
use sort_module
!
    implicit none
!
#include "asterc/asmpi_comm.h"
#include "asterf_types.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "jeveux.h"
!
    integer, intent(in) :: domdist(*), nbdom
    integer, intent(out) :: comm(*), tag(*)
!
!---------------------------------------------------------------------------------------------------
!
! Le but est de construire un arbre de comm optimisé pour les comm point à point
!
!---------------------------------------------------------------------------------------------------
!
#ifdef ASTER_HAVE_MPI
!
    integer :: rank, nbproc, nbdist_tot, i_proc, max_nbdom, nb_comm, nbdom_inf, domtmp(50)
    integer :: dom1, dom2, i_dom, ind, nb_comm_loc
    mpi_int :: mrank, msize, mpicou, count_send, count_recv
    aster_logical :: find
    integer, pointer :: v_nbdist(:) => null()
    integer, pointer :: v_deca(:) => null()
    integer, pointer :: v_dist(:) => null()
    integer, pointer :: v_send(:) => null()
    integer, pointer :: v_recv(:) => null()
    integer, pointer :: v_rest(:) => null()
    aster_logical, pointer :: v_proc(:) => null()

!
    call jemarq()
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank = mrank, size = msize)
    rank = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    ASSERT(nbdom < 100)
    nb_comm = 0
!
    nbdom_inf = 0
    domtmp = -1
    do i_dom = 1, nbdom
        if( domdist(i_dom) > rank) then
            nbdom_inf = nbdom_inf + 1
            domtmp(nbdom_inf) = domdist(i_dom)
        end if
    end do
    call sort_i8(domtmp, nbdom_inf)
! --- On compte le nombre de sous-domaine
    AS_ALLOCATE(vi=v_nbdist, size=nbproc)
    count_send = to_mpi_int(1)
    count_recv = count_send
    call asmpi_allgather_i([nbdom_inf], count_send, v_nbdist, count_recv, mpicou)
    max_nbdom = maxval(v_nbdist)
    AS_ALLOCATE(vi=v_deca, size=nbproc)
!
    v_deca = 0
    do i_proc = 1, nbproc-1
        v_deca(i_proc+1) = v_deca(i_proc) + v_nbdist(i_proc)
    end do
    nbdist_tot = v_deca(nbproc) + v_nbdist(nbproc)
    if( nbdist_tot == 0 .or. max_nbdom == 0) go to 999
!
! --- On récupère la liste des sous-domaines
    AS_ALLOCATE(vi=v_send, size=max_nbdom)
    AS_ALLOCATE(vi=v_recv, size=nbproc*max_nbdom)
    v_send = 0
    v_send(1:nbdom_inf) = domtmp(1:nbdom_inf)
    count_send = to_mpi_int(max_nbdom)
    count_recv = count_send
    call asmpi_allgather_i(v_send, count_send, v_recv, count_recv , mpicou)
    AS_ALLOCATE(vi=v_dist, size=nbdist_tot)
!
    do i_proc = 1, nbproc
        v_dist(v_deca(i_proc)+1:v_deca(i_proc+1)) = &
            v_recv((i_proc-1)*max_nbdom+1:(i_proc-1)*max_nbdom+v_nbdist(i_proc))
    end do
!
! --- On crée la liste des comm
    AS_ALLOCATE(vl=v_proc, size=nbproc)
    AS_ALLOCATE(vi=v_rest, size=nbproc)
    v_rest = v_nbdist
    nb_comm = 0
    nb_comm_loc = 0
    do while( nb_comm < nbdist_tot)
        v_proc = ASTER_TRUE
        do i_proc = 1, nbproc
            find = ASTER_FALSE
            if(v_proc(i_proc) .and. v_rest(i_proc) > 0) then
                dom1 = i_proc - 1
                do i_dom = 1, v_nbdist(i_proc)
                    ind = v_deca(i_proc)+i_dom
                    dom2 = v_dist(ind)
                    if(dom2 .ne. -1 .and. v_proc(dom2)) then
                        nb_comm = nb_comm + 1
                        v_proc(dom1+1) = ASTER_FALSE
                        v_proc(dom2+1) = ASTER_FALSE
                        v_rest(i_proc) = v_rest(i_proc) - 1
                        v_dist(ind) = -1
                        if(dom1 == rank) then
                            nb_comm_loc = nb_comm_loc + 1
                            tag(nb_comm_loc) = nb_comm
                            comm(nb_comm_loc) = dom2
                        elseif(dom2 == rank) then
                            nb_comm_loc = nb_comm_loc + 1
                            tag(nb_comm_loc) = nb_comm
                            comm(nb_comm_loc) = dom1
                        end if
                        !if(rank == 0) print*, dom1, dom2, v_rest(i_proc)
                        exit
                    end if
                end do
            end if
        end do
    end do
!
    ASSERT(nb_comm_loc == nbdom)
!
    AS_DEALLOCATE(vl=v_proc)
    AS_DEALLOCATE(vi=v_rest)
    AS_DEALLOCATE(vi=v_send)
    AS_DEALLOCATE(vi=v_recv)
    AS_DEALLOCATE(vi=v_dist)
999 continue
    AS_DEALLOCATE(vi=v_deca)
    AS_DEALLOCATE(vi=v_nbdist)
    call jedema()
#endif
!
end subroutine
