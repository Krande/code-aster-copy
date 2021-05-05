! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine crnlgc(numddl)
    implicit none
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_i.h"
#include "asterc/asmpi_send_i.h"
#include "asterc/loisem.h"
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/utimsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=14) :: numddl

#ifdef ASTER_HAVE_MPI
!
    integer :: rang, nbproc, jrefn, nbjoin, iaux, nddll
    integer :: iproc1, iproc2, posit, nmatch, icmp, ico2, nbcmp
    integer :: dime, idprn1, idprn2, ntot, num, ili, nunoel
    integer :: nec, l, numpro, jjoine, jjoinr, nbnoee, jaux, numno1, numno2, iec
    integer :: ncmpmx, iad, jcpnec, jencod, jenvoi1, lgenve1, lgenvr1, poscom, nddld
    integer :: nbddll, jnequ, nddl, jenco2, deccmp, jcpne2
    integer :: jnbddl, decalp, jddlco, posddl, nuddl, inttmp, nbnoer, jrecep1
    integer :: decalm, nbjver, jprddl, jnujoi, cmpteu, lgrcor, jnbjoi, curpos
    integer :: jmlogl, nuno, ieq, numero_noeud, nb_ddl_envoi, nbddl
    integer :: ibid, nddlg, jenvoi2, jrecep2, jjoin, ijoin, numnoe
    integer :: ifm, niv, vali(5), ino, nno, nb_node, nlag
    integer :: lgenve2, lgenvr2, jnujoi1, jnujoi2, iret, iret1, iret2, nlili
    integer(kind=4) :: un
    real(kind=8) :: dx, dy, dz
    integer(kind=4) :: iaux4, num4, numpr4, n4e, n4r
    mpi_int :: mrank, mnbproc, mpicou
    integer, pointer :: v_noex(:) => null()
    integer, pointer :: v_nugll(:) => null()
    integer, pointer :: v_posdd(:) => null()
    integer, pointer :: v_dojoin(:) => null()
    integer, pointer :: v_graloc(:) => null()
    integer, pointer :: v_gracom(:) => null()
    integer, pointer :: v_mask(:) => null()
    integer, pointer :: v_tmp(:) => null()
    integer, pointer :: v_ordjoi(:) => null()
    integer, pointer :: v_nbjo(:) => null()
    integer, pointer :: v_deeq(:) => null()
!
    character(len=4) :: chnbjo
    character(len=8) :: mesh, k8bid, nomgdr
    character(len=19) :: nomlig
    character(len=24) :: nojoie, nojoir, nomtmp, mesh24, nonulg, join
!
!----------------------------------------------------------------------
    integer :: zzprno
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
    zzprno(ili,nunoel,l) = zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
    call jemarq()
!
    call infniv(ifm,niv)
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank = mrank, size = mnbproc)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(mnbproc)
    ASSERT(nbproc.lt.9999)

!   RECUPERATION DU NOM DU MAILLAGE DANS LE BUT D'OBTENIR LE JOINT
    call jeveuo(numddl//'.NUME.REFN', 'L', jrefn)
    mesh = zk24(jrefn)
    nomgdr = zk24(jrefn + 1)

    call jeveuo(numddl//'.NUME.NULG', 'E', vi=v_nugll)
    call jeveuo(numddl//'.NUME.PDDL', 'E', vi=v_posdd)

    call wkvect('&&CRNULG.GRAPH_COMM', 'V V I', nbproc*nbproc, vi=v_gracom)
    call wkvect('&&CRNULG.GRAPH_LOC', 'V V I', nbproc, vi=v_graloc)
!
! -- On crée la liste des proc avec lequel rank doit communiquer - 1 si oui
    nbjoin = 0
    call jeexin(mesh//'.DOMJOINTS', iret)
    if(iret > 0) then
        call jeveuo(mesh//'.DOMJOINTS', 'L', vi=v_dojoin)
        call jelira(mesh//'.DOMJOINTS', 'LONMAX', nbjoin, k8bid)
    !   CREATION DU GRAPH LOCAL
        do iaux = 1, nbjoin
            v_graloc(v_dojoin(iaux)+1) = 1
        end do
        nbjoin = nbjoin/2
    endif

!   ON COMMUNIQUE POUR SAVOIR QUI EST EN RELATION AVEC QUI
    call asmpi_allgather_i(v_graloc, mnbproc, v_gracom, mnbproc, mpicou)

!   RECHERCHE DES COUPLAGES MAXIMAUX
    call wkvect('&&CRNULG.MASQUE', 'V V I', nbproc*nbproc, vi=v_mask)
    nmatch = 1

    do iproc1 = 1, nbproc
        do iproc2 = 1, nbproc
            posit = (iproc1-1)*nbproc + iproc2
            if (v_gracom(posit) == 1) then
                v_gracom(posit) = 0
                v_mask(posit) = nmatch
                posit = (iproc2-1)*nbproc + iproc1
                v_gracom(posit) = 0
                v_mask(posit) = nmatch
                nmatch = nmatch + 1
            endif
        end do
    end do

    call jedetr('&&CRNULG.GRAPH_COMM')
    call jedetr('&&CRNULG.GRAPH_LOC')

    nmatch = nmatch - 1
    call wkvect('&&CRNULG.ORDJOI', 'V V I', nmatch, vi=v_ordjoi)
    v_ordjoi(:) = -1

    nbjver = 0
    do iaux = 0, nbproc - 1
        num = v_mask(rang*nbproc + iaux + 1)
        ASSERT(num .le. nmatch)
        if (num .ne. 0) then
            v_ordjoi(num) = iaux
            nbjver = nbjver + 1
        endif
    end do
    ASSERT(nbjver .eq. nbjoin)
    call wkvect(numddl//'.NUME.NBJO', 'G V I', 1 + nmatch, vi=v_nbjo)
    v_nbjo(1) = nmatch

!     !!!! IL PEUT ETRE INTERESSANT DE STOCKER CES INFOS
!     !!!! EN CAS DE CONSTRUCTION MULTIPLE DE NUMEDDL
!
!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(numddl//'.NUME.PRNO', 'E', idprn1)
    call jeveuo(jexatr(numddl//'.NUME.PRNO', 'LONCUM'), 'L', idprn2)
    call jelira(jexnum(numddl//'.NUME.PRNO', 1), 'LONMAX', ntot, k8bid)
    call jeveuo(numddl//'.NUME.DEEQ', 'L', vi=v_deeq)

    call jeveuo(mesh//'.DIME', 'L', dime)

!   !!! VERIFIER QU'IL N'Y A PAS DE MACRO-ELTS
!   CALCUL DU NOMBRE D'ENTIERS CODES A PARTIR DE LONMAX
    nec = ntot/zi(dime) - 2
    call wkvect('&&CRNULG.NEC', 'V V I', nec, jencod)
    call wkvect('&&CRNULG.NEC2', 'V V I', nec, jenco2)

    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgdr), 'L', iad)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgdr), 'LONMAX', ncmpmx, k8bid)
    call wkvect('&&CRNULG.CMP', 'V V I', ncmpmx, jcpnec)
    call wkvect('&&CRNULG.CMP2', 'V V I', ncmpmx, jcpne2)
!
!   Il faut maintenant communiquer les numeros partages
!   NOTE : On pourrait sans doute se passer d'une communication puisque celui qui recoit
!          sait ce qu'il attend et celui qui envoit pourrait tout envoyer
    do iaux = 0, nmatch - 1
!
!       RECHERCHE DU JOINT
        numpro = v_ordjoi(iaux+1)
        v_nbjo(iaux + 2) = numpro
        if (numpro .ne. -1) then
        num = v_mask(rang*nbproc + numpro + 1)
        call codent(numpro, 'G', chnbjo)
        nojoie = mesh//'.E'//chnbjo
        nojoir = mesh//'.R'//chnbjo
        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
        call jeveuo(nojoir, 'L', jjoinr)
        call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
        nbnoee = nbnoee/2
        nbnoer = nbnoer/2
!
!       DES DEUX COTES LES NOEUDS NE SONT PAS DANS LE MEME ORDRE ?
        num4 = num
        numpr4 = numpro
        lgenve1 = nbnoee*(1 + nec) + 1
        lgenvr1 = nbnoer*(1 + nec) + 1
        call wkvect('&&CRNULG.NOEUD_NEC_E1', 'V V I', lgenvr1, jenvoi1)
        call wkvect('&&CRNULG.NOEUD_NEC_R1', 'V V I', lgenve1, jrecep1)
!
        lgenve2 = nbnoee*(1 + nec) + 1
        lgenvr2 = nbnoer*(1 + nec) + 1
!
!       On commence par envoyer, le but final est de recevoir les numeros de ddl
!       On boucle donc sur les noeuds a recevoir
        nb_ddl_envoi = 0
        do jaux = 1, nbnoer
            poscom = (jaux - 1)*(1 + nec) + 1
            numno1 = zi(jjoinr + 2*(jaux - 1))
            numno2 = zi(jjoinr + 2*jaux - 1)
            zi(jenvoi1 + poscom) = numno2
            do iec = 1, nec
                zi(jenvoi1 + poscom + iec) = zzprno(1, numno1, 2 + iec)
            end do
            nb_ddl_envoi = nb_ddl_envoi + zzprno(1, numno1, 2)
        end do
        zi(jenvoi1) = nb_ddl_envoi
        n4r = lgenvr1
        n4e = lgenve1
        if (rang .lt. numpro) then
            call asmpi_send_i(zi(jenvoi1), n4r, numpr4, num4, mpicou)
            call asmpi_recv_i(zi(jrecep1), n4e, numpr4, num4, mpicou)
        else if (rang.gt.numpro) then
            call asmpi_recv_i(zi(jrecep1), n4e, numpr4, num4, mpicou)
            call asmpi_send_i(zi(jenvoi1), n4r, numpr4, num4, mpicou)
        else
            ASSERT(.false.)
        endif
        call wkvect('&&CRNULG.NUM_DDL_GLOB_E', 'V V I', zi(jrecep1), jenvoi2)
        call wkvect('&&CRNULG.NUM_DDL_GLOB_R', 'V V I', zi(jenvoi1), jrecep2)
        call codent(iaux, 'G', chnbjo)
        call wkvect(numddl//'.NUMEE'//chnbjo, 'G V I', zi(jrecep1), jnujoi1)
!
        nbddl = 0
        do jaux = 1, nbnoee
            poscom = (jaux - 1)*(1 + nec) + 1
            numno1 = zi(jrecep1 + poscom)
!
            nddl = zzprno(1, numno1, 1)
            nddlg = v_nugll(nddl)
!
!           Recherche des composantes demandees
            do iec = 1, nec
                zi(jencod + iec - 1) = zzprno(1, numno1, 2 + iec)
                zi(jenco2 + iec - 1) = zi(jrecep1 + poscom + iec)
            end do
            call isdeco(zi(jencod), zi(jcpnec), ncmpmx)
            call isdeco(zi(jenco2), zi(jcpne2), ncmpmx)
            ico2 = 0
            do icmp = 1, ncmpmx
                if ( zi(jcpnec + icmp - 1) .eq. 1 ) then
                    if ( zi(jcpne2 + icmp - 1) .eq. 1 ) then
                        ASSERT(nddlg .ne. -1)
                        zi(jenvoi2 + nbddl) = nddlg + ico2
                        zi(jnujoi1 + nbddl) = nddl + ico2
                    endif
                    ico2 = ico2 + 1
                    nbddl = nbddl + 1
                endif
            enddo
        enddo
!
        ASSERT(zi(jrecep1) .eq. nbddl)
        n4e = nbddl
        n4r = nb_ddl_envoi
        if (rang .lt. numpro) then
            call asmpi_send_i(zi(jenvoi2), n4e, numpr4, num4, mpicou)
            call asmpi_recv_i(zi(jrecep2), n4r, numpr4, num4, mpicou)
        else if (rang.gt.numpro) then
            call asmpi_recv_i(zi(jrecep2), n4r, numpr4, num4, mpicou)
            call asmpi_send_i(zi(jenvoi2), n4e, numpr4, num4, mpicou)
        else
            ASSERT(.false.)
        endif
        call wkvect(numddl//'.NUMER'//chnbjo, 'G V I', nb_ddl_envoi, jnujoi2)
!
        curpos = 0
        do jaux = 1, nbnoer
            numno1 = zi(jjoinr + 2*(jaux - 1))
            nddll = zzprno(1, numno1, 1)
            nbcmp = zzprno(1, numno1, 2)
            do icmp = 0, nbcmp - 1
                ASSERT(zi(jrecep2 + curpos) .ne. -1)
                v_nugll(nddll + icmp) = zi(jrecep2 + curpos)
                v_posdd(nddll + icmp) = numpro
                zi(jnujoi2 + curpos) = nddll + icmp
                curpos = curpos + 1
            enddo
        enddo
        ASSERT(curpos .eq. nb_ddl_envoi)
!
        call jedetr('&&CRNULG.NOEUD_NEC_E1')
        call jedetr('&&CRNULG.NOEUD_NEC_R1')
        call jedetr('&&CRNULG.NUM_DDL_GLOB_E')
        call jedetr('&&CRNULG.NUM_DDL_GLOB_R')
        endif
    end do

    call jelira(numddl//'.NUME.PRNO', 'NMAXOC', nlili, k8bid)
    do ili = 2, nlili
        call jeexin(jexnum(numddl//'.NUME.PRNO', ili), iret)
        if( iret.ne.0 ) then
            call jenuno(jexnum(numddl//'.NUME.LILI', ili), nomlig)
            join = nomlig//".NBJO"
            call jeexin(join, iret)
            if( iret.eq.0 ) cycle
            call jeveuo(join, 'L', jjoin)
            call jelira(join, 'LONMAX', nbjoin)
            do ijoin = 1, nbjoin
                numpro = zi(jjoin+ijoin-1)
                if( numpro.ne.-1 ) then
                    numpr4 = numpro
                    num4 = ijoin
                    call codent(numpro, 'G', chnbjo)
                    nojoie = nomlig//'.E'//chnbjo
                    nojoir = nomlig//'.R'//chnbjo

                    call jeexin(nojoie, iret1)
                    if( iret1.ne.0 ) then
                        call jeveuo(nojoie, 'L', jjoine)
                        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_E', 'V V I', nbnoee, jnujoi1)
                        do jaux = 1, nbnoee
                            numnoe = -zi(jjoine+jaux-1)
                            ASSERT(zzprno(ili, numnoe, 2).eq.1)
                            nddll = zzprno(ili, numnoe, 1)
                            zi(jnujoi1+jaux-1) = v_nugll(nddll)
                        enddo
                        n4e = nbnoee
                    endif

                    call jeexin(nojoir, iret2)
                    if( iret2.ne.0 ) then
                        call jeveuo(nojoir, 'L', jjoinr)
                        call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_R', 'V V I', nbnoer, jnujoi2)
                        n4r = nbnoer
                    endif

                    if (rang .lt. numpro) then
                        if( iret1.ne.0 ) then
                            call asmpi_send_i(zi(jnujoi1), n4e, numpr4, num4, mpicou)
                        endif
                        if( iret2.ne.0 ) then
                            call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, num4, mpicou)
                        endif
                    else if (rang.gt.numpro) then
                        if( iret2.ne.0 ) then
                            call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, num4, mpicou)
                        endif
                        if( iret1.ne.0 ) then
                            call asmpi_send_i(zi(jnujoi1), n4e, numpr4, num4, mpicou)
                        endif
                    endif

                    if( iret2.ne.0 ) then
                        do jaux = 1, nbnoer
                            numnoe = -zi(jjoinr+jaux-1)
                            ASSERT(zzprno(ili, numnoe, 2).eq.1)
                            nddll = zzprno(ili, numnoe, 1)
                            v_nugll(nddll) = zi(jnujoi2+jaux-1)
                            v_posdd(nddll) = numpro
                        enddo
                    endif
                    call jedetr('&&CRNUGL.NUM_DDL_GLOB_E')
                    call jedetr('&&CRNUGL.NUM_DDL_GLOB_R')
                endif
            enddo
        endif
    enddo
!
    call jedetc('V', '&&CRNULG', 1)
!
! --- Vérification de la numérotation
!
!   NOMBRE DE DDL LOCAUX
    call jeveuo(numddl//'.NUME.NEQU', 'L', jnequ)
    nbddll = zi(jnequ)
    call jeveuo(mesh//'.NOEX', 'L', vi=v_noex)
    do iaux = 1, nbddll
        nuno = v_deeq((iaux-1)*2+1)
        if(nuno.ne.0) then
            ASSERT(v_posdd(iaux) == v_noex(nuno))
        else
            ASSERT(v_posdd(iaux) .ne. -1)
        end if
    end do
!
! --- Affichage inconnue systeme
    if (niv .ge. 1 ) then
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node)
        nno = 0
        do ino = 1, nb_node
            if(v_noex(ino) == rang) then
                if (zi(idprn1-1+ (ino-1)* (2+nec)+2) .gt. 0) nno = nno + 1
            end if
        end do
        call asmpi_comm_vect("MPI_SUM", "I", sci=nno)
!
        call jeexin(numddl//'.NUME.MDLA', iret)
        if(iret .ne. 0) then
            call jelira(numddl//'.NUME.MDLA', 'LONMAX', nlag, k8bid)
            nlag = nlag / 3
        else
            nlag = 0
        end if
        call asmpi_comm_vect("MPI_SUM", "I", sci=nlag)
!
        vali(1) = zi(jnequ + 1)
        vali(2) = zi(jnequ + 1) - nlag
        vali(3) = nno
        vali(4) = nlag
        vali(5) = nlag / 2
        call utmess('I', 'FACTOR_1', ni=5, vali=vali)
    end if
!
! --- Pour debuggage en hpc
    if(ASTER_FALSE) then
        nonulg = mesh//'.NULOGL'
        call jeveuo(nonulg, 'L', jmlogl)
        do iaux = 0, nbddll - 1
            nuno = v_deeq(iaux*2+1)
            if(nuno.ne.0) nuno = zi(jmlogl + nuno - 1) + 1
! numero ddl local, numéro noeud local, numéro noeud global, num composante du noeud,
!            num ddl global, num proc proprio
            write(130+rang, *) iaux, v_deeq(iaux*2+1), nuno , v_deeq(iaux*2 + 2), &
            v_nugll(iaux + 1), v_posdd(1 + iaux)
        end do
        flush(130+rang)
    end if
!
    call jedetr('&&CRNULG.ORDJOI')
!
    call jedema()
#else
    character(len=14) :: kbid
    kbid = numddl
#endif
!
end subroutine
