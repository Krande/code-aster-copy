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
subroutine lrmjoi(fid, nommail, nomam2, nbnoeu)
!
use sort_module
!
    implicit none
#include "asterf.h"
#include "asterf_med.h"
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterc/ismaem.h"
#include "asterfort/as_mmhgnr.h"
#include "asterfort/as_msdcrr.h"
#include "asterfort/as_msdjni.h"
#include "asterfort/as_msdnjn.h"
#include "asterfort/as_msdszi.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/decode_join.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lrm_clean_joint.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
    med_idt, intent(in) :: fid
    integer, intent(in) :: nbnoeu
    character(len=*), intent(in) :: nomam2, nommail
!
! ---------------------------------------------------------------------------------------------
!
! LECTURE DU FORMAT MED : Lecture des joints pour le ParallelMesh
!
! ---------------------------------------------------------------------------------------------
!
    character(len=4) :: chdomdis
    character(len=7) :: code
    character(len=8) :: mesh, nom
    character(len=24) :: nonulg, nojoin, connex, nomnoe
    character(len=MED_NAME_SIZE) :: nomjoi, nommad
    character(len=MED_COMMENT_SIZE) :: descri
    integer :: rang, nbproc, nbjoin, domdis, nstep, ncorre
    integer :: icor, entlcl, geolcl, entdst, geodst, ncorr2
    integer :: jnlogl, codret, i_join, ino, numno, deca
    integer :: ima, nbma, node_id, nbnoma, dom1, dom2, incr, nbjoin_max
    mpi_int :: mrank, msize
    integer, pointer :: v_noext(:) => null()
    integer, pointer :: v_nojoin(:) => null()
    integer, pointer :: v_maex(:) => null()
    integer, pointer :: v_connex(:) => null()
    integer, pointer :: v_dom(:) => null()
!
    call jemarq()
!
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    ASSERT(nbproc <= MT_DOMMAX)
!
! --- Uniquement pour les ParallelMesh
!
    mesh = nommail(1:8)
!
! --- L'objet .NULOGL permet d'avoir la num??rotation globale des noeuds
!
    nonulg = mesh//'.NULOGL'
    call wkvect(nonulg, 'G V I', nbnoeu, jnlogl)
!
! --- On renomme les noeuds avec le num??ro global
!
    nomnoe = mesh//'.NOMNOE'
    call jedetr(nomnoe)
    call jecreo(nomnoe, 'G N K8')
    call jeecra(nomnoe, 'NOMMAX', nbnoeu)
    do ino = 1 , nbnoeu
        call codlet(ino, 'G', code)
        nom = 'N'//code
        call jecroc(jexnom(nomnoe, nom))
    end do
!
! --- L'objet .NOEX permet de savoir ?? qui appartient le noeud.
!     Pour les noeuds internes, c'est le proc courant
!     Pour les noeuds joints, c'est un autre proc et il faut le trouver par lecture des joints
!     Par d??faut, tout les noeuds d'un domaine appartient au moins ?? ce domaine
!
    call wkvect(mesh//'.NOEX', 'G V I', nbnoeu, vi=v_noext)
    v_noext(1:nbnoeu) = rang
!
    if ( nbproc > 1 ) then
!
! --- R??cup??ration de la num??rotation globale des noeuds
!
        call as_mmhgnr(fid, nomam2, MED_NODE, MED_NONE, zi(jnlogl), nbnoeu, codret)
!
! --- R??cup??ration du nombre de joints
!
        call as_msdnjn(fid, nomam2, nbjoin, codret)
!
! --- On lit l'info de tout les joints quelques soient leurs types car med le permet
!     mais il faudra faire un tri apr??s pour garder que ceux qui nous int??resse
!     Type de joints (entre deux noeuds uniquement pour le moment):
!     - noeud interne - noeud joint (il faut le garder)
!     - noeud joint   - noeud joint (il faut le supprimer car ne sert ?? rien pour nous)
!
!     Le probl??me c'est que l'on ne sais pas ?? l'avance qui est un noeud interne
!     et qui est un noeud joint. C'est l'objet .NOEX qui l'indique mais il faut faire
!     des comm pour le savoir
!
        nbjoin_max = nbjoin
        call asmpi_comm_vect("MPI_MAX", "I", sci=nbjoin_max)
        if(nbjoin_max == 0) then
            call utmess('A', 'MAILLAGE1_5', sk=nomam2)
        end if

        if(nbjoin > 0) then
            call wkvect(mesh//'.DOMJOINTS', 'G V I', nbjoin/2, vi=v_dom)
!
! --- Boucle sur les joints entre les sous-domaines
!
            incr = 0
            do i_join = 1, nbjoin
                call as_msdjni(fid, nomam2, i_join, nomjoi, descri, domdis, &
                            nommad, nstep, ncorre, codret)
                ASSERT(domdis <= nbproc)
                ASSERT(ncorre == 1)
!
                do icor = 1, ncorre
                    call as_msdszi(fid, nomam2, nomjoi, MED_NO_DT, MED_NO_IT, icor, entlcl, &
                                geolcl, entdst, geodst, ncorr2, codret)

                    call decode_join(nomjoi, dom1, dom2)
!
                    if ( entlcl.eq.MED_NODE.and.geolcl.eq.MED_NONE ) then
                        call codlet(domdis, 'G', chdomdis)
                        if ( dom1.eq.rang ) then
                            nojoin = mesh//'.R.'//chdomdis
                            incr = incr + 1
                            v_dom(incr) = domdis
                        else
                            nojoin = mesh//'.E.'//chdomdis
                        endif
!
! --- R??cup??ration de la table de correspondance pour les noeuds partag??s par 2 sous-domaines
!
                        call wkvect(nojoin, 'V V I', 2*ncorr2, vi=v_nojoin)
                        call as_msdcrr(fid, nomam2, nomjoi, -1, -1, entlcl, &
                                    geolcl, entdst, geodst, 2*ncorr2, &
                                    v_nojoin, codret)
!
! --- On r??cup??re le num??ro du sous-domaine pour les noeuds partag??s
!
                        if(dom1.eq.rang) then
                            deca = 1
                            do ino = 1, ncorr2
                                numno = v_nojoin(deca)
                                v_noext(numno) = -1
                                deca = deca +2
                            end do
                        end if
!
                    endif
                enddo
            enddo
!
            call sort_i8(v_dom, nbjoin/2)
        end if
!
! --- On nettoie les joints des noeuds en trop
!
        call lrm_clean_joint(mesh, v_noext)
!
! --- Verification NOEX
!
        do ino = 1, nbnoeu
            ASSERT(v_noext(ino).ne.-1)
        end do
!
    else
        do ino = 1, nbnoeu
            zi(jnlogl + ino - 1) = ino - 1
        end do
    endif
!
! --- Creation .MAEX
!
    connex = mesh //'.CONNEX'
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
    call wkvect(mesh//'.MAEX', 'G V I', nbma, vi=v_maex)
    v_maex(1:nbma) = ismaem()
    do ima = 1, nbma
        call jelira(jexnum(connex , ima), 'LONMAX', nbnoma)
        call jeveuo(jexnum(connex , ima), 'L', vi=v_connex)
        do ino = 1, nbnoma
            node_id = v_connex(ino)
            v_maex(ima) = min(v_maex(ima), v_noext(node_id))
        end do
    end do
!
    call jedema()
!
end subroutine
