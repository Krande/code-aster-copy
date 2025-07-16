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
!
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine lrm_clean_joint(mesh, v_noex)
!
    implicit none
#include "asterf.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecrec.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
    character(len=8) :: mesh
    integer(kind=8), intent(inout) :: v_noex(*)
!
!
! ---------------------------------------------------------------------------------------------
!
! LECTURE DU FORMAT MED : Nettoyage des joints pour le ParallelMesh
!
!     On a lu l'info de tout les joints quelques soient leurs types car med le permet
!     mais il faut faire un tri maintenant dans les correspondances pour garder que celles
!     qui nous intéresse
!     Type de correspondance (entre deux noeuds uniquement pour le moment):
!     - noeud interne - noeud joint (il faut le garder)
!     - noeud joint   - noeud joint (il faut le supprimer car ne sert à rien pour nous)
!
! ---------------------------------------------------------------------------------------------
!
    character(len=4) :: chdomdis
    character(len=8) :: k8bid
    character(len=19) :: comm_name, tag_name, joints
    character(len=24) :: name_join_e_old, send, name_join_r_old, recv, gcom, pgid
    integer(kind=8) :: rang, domdis, nbproc, i_comm, nb_comm, domj_i, domdi2
    integer(kind=8) :: nb_corr, ino, numno, deca
    integer(kind=8) :: nb_node_e, nb_node_r
    mpi_int :: mpicou, count_send, count_recv, tag, id, mrank, msize
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_nojoe(:) => null()
    integer(kind=8), pointer :: v_nojor(:) => null()
    integer(kind=8), pointer :: v_domj(:) => null()
    integer(kind=8), pointer :: v_gcom(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
    integer(kind=8), pointer :: v_name_join_e_old(:) => null()
    integer(kind=8), pointer :: v_name_join_e_new(:) => null()
    integer(kind=8), pointer :: v_name_join_r_old(:) => null()
    integer(kind=8), pointer :: v_name_join_r_new(:) => null()
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank=mrank, size=msize)
    if (msize .le. 1) then
!       TODO: rename ".E." to ".E"...
        goto 999
    end if
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    DEBUG_MPI('lrm_clean_joint', rang, nbproc)
!
    call jemarq()
!
! --- Create COMM_GRAPH
    comm_name = '&&LRMCLEAN.COMM'
    tag_name = '&&LRMCLEAN.TAG'
    call create_graph_comm(mesh, "MAILLAGE_P", nb_comm, comm_name, tag_name)

!   -- creation de la collection dipersee .JOIN  :
    joints = mesh//".JOIN"
    send = joints//".SEND"
    recv = joints//".RECV"
    gcom = joints//".GCOM"
    pgid = joints//".PGID"
    if (nb_comm > 0) then
        call jecrec(send, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_comm)
        call jecrec(recv, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_comm)
        call jeveuo(joints//".DOMJ", 'L', vi=v_domj)
        call jeveuo(gcom, 'L', vi=v_gcom)
        call jeveuo(pgid, 'L', vi4=v_pgid)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
        mpicou = v_gcom(1)
    end if

    do i_comm = 1, nb_comm
        domj_i = v_comm(i_comm)
        domdis = v_domj(domj_i)
        domdi2 = v_pgid(domdis+1)
        call codlet(domdis, 'G', chdomdis)
!
! --- Il faut préparer les noeuds à envoyer et à recevoir
!
        name_join_e_old = mesh//".E."//chdomdis
        call jelira(name_join_e_old, 'LONMAX', nb_corr, k8bid)
        nb_node_e = nb_corr/2
!
        call wkvect("&&LRMJOI.NOJOE", 'V V I', nb_node_e, vi=v_nojoe)

        name_join_r_old = mesh//".R."//chdomdis
        call jelira(name_join_r_old, 'LONMAX', nb_corr, k8bid)
        nb_node_r = nb_corr/2
!
        call wkvect("&&LRMJOI.NOJOR", 'V V I', nb_node_r, vi=v_nojor)
!
! --- 0 le noeud m'appartient et -1 le noeud ne m'appartient pas de .ET
!
        call jeveuo(name_join_e_old, 'L', vi=v_name_join_e_old)
        deca = 1
        nb_corr = 0
        do ino = 1, nb_node_e
            numno = v_name_join_e_old(deca)
            if (v_noex(numno) .ne. rang) then
                v_nojoe(ino) = -1
            else
                nb_corr = nb_corr+1
            end if
            deca = deca+2
        end do
        ASSERT(nb_corr > 0)
!
! --- On crée le nouveau joint .E
!
        call jecroc(jexnum(send, i_comm))
        call jeecra(jexnum(send, domj_i), 'LONMAX', 2*nb_corr)
        call jeecra(jexnum(send, domj_i), 'LONUTI', 2*nb_corr)
        call jeveuo(jexnum(send, domj_i), 'E', vi=v_name_join_e_new)
!
        deca = 0
        do ino = 1, nb_node_e
            if (v_nojoe(ino) == 0) then
                v_name_join_e_new(deca+1) = v_name_join_e_old((ino-1)*2+1)
                v_name_join_e_new(deca+2) = v_name_join_e_old((ino-1)*2+2)
                deca = deca+2
            end if
        end do
        ASSERT(deca == 2*nb_corr)
!
! --- Envoie des informations pour le joint
!
        tag = to_mpi_int(v_tag(i_comm))
        id = to_mpi_int(domdi2)
        count_send = to_mpi_int(nb_node_e)
        count_recv = to_mpi_int(nb_node_r)
        call asmpi_sendrecv_i(v_nojoe, count_send, id, tag, &
                              v_nojor, count_recv, id, tag, mpicou)
!
! --- Il faut supprimer les correspondances en trop maintenant dans le joint .RT
!
        nb_corr = 0
        do ino = 1, nb_node_r
            if (v_nojor(ino) == 0) then
                nb_corr = nb_corr+1
            end if
        end do
        ASSERT(nb_corr > 0)
!
        call jecroc(jexnum(recv, i_comm))
        call jeecra(jexnum(recv, domj_i), 'LONMAX', 2*nb_corr)
        call jeecra(jexnum(recv, domj_i), 'LONUTI', 2*nb_corr)
        call jeveuo(jexnum(recv, domj_i), 'E', vi=v_name_join_r_new)
        call jeveuo(name_join_r_old, 'L', vi=v_name_join_r_old)
!
        deca = 0
        do ino = 1, nb_node_r
            if (v_nojor(ino) == 0) then
                v_name_join_r_new(deca+1) = v_name_join_r_old((ino-1)*2+1)
                v_name_join_r_new(deca+2) = v_name_join_r_old((ino-1)*2+2)
                ASSERT(v_noex(v_name_join_r_new(deca+1)) == -1)
                v_noex(v_name_join_r_new(deca+1)) = domdis
                deca = deca+2
            end if
        end do
        ASSERT(deca == 2*nb_corr)
!
        call jedetr("&&LRMJOI.NOJOE")
        call jedetr("&&LRMJOI.NOJOR")
        call jedetr(name_join_e_old)
        call jedetr(name_join_r_old)
    end do
!
! --- Il faut détruire les objets temporaires
!
    call jedetr(comm_name)
    call jedetr(tag_name)
!
    call jedema()
999 continue

end subroutine
