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

subroutine crnustd(numddl)
    implicit none
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_i.h"
#include "asterc/asmpi_send_i.h"
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=14) :: numddl
!
! --------------------------------------------------------------------------------------------------
!
!  Le but est de créer une numérotation parallèle qui soit toujours la même quelque soit le
!  nombre de sous-domaines
!
!  .NUME.DEEG(2*nbbdl) : 2*(i_ddl-1) + 1 -> numéro global du noeud physique ou tardif
!                      : 2*(i_ddl-1) + 2 -> numéro de la composante (comme DEEQ)
!
!  .NUME.NULS(nbddl) : iddl -> numéro d'équation globale (différent de .NUME.NULG)
!    qui ne dépend pas de la partition du maillage
!
! --------------------------------------------------------------------------------------------------
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
!
    integer :: ili, nunoel, idprn1, idprn2, ntot, lonmax, nbno_prno
    integer :: nbddll, i_proc, ino, iret, nbcmp, iec, iret1, iret2, jjoin, jjoine
    integer :: numero_noeud, numero_cmp, rang, nbproc, jrefn
    integer :: nec, numloc, dime, nbddl_lag, i_ddl, nddl, nddlg, nddll
    integer :: nbno, nbno_lc, nbno_gl, nbno_max, nbddll_gl, numnoe
    integer :: nbddl_phys_gl, nbddl_lag_gl, iproc, nbjoin, i_join, jnujoi1, jnujoi2
    integer :: numpro, nbnoee, nbnoer, jjoinr, poscom, numno1, numno2
    integer :: lgenve1, lgenvr1, jencod, jenco2, lgenve2, lgenvr2
    integer :: jaux, nb_ddl_envoi, jrecep1, jenvoi1, jenvoi2, jrecep2
    integer :: nbddl, ncmpmx, iad, jcpnec, jcpne2, ico2, icmp, curpos
    integer :: nbno_lili_lc, nbno_lili_gl
    mpi_int :: mrank, msize, mpicou, nbno4
    mpi_int :: iaux4, num4, numpr4, n4e, n4r
    mpi_int, parameter :: one4 = to_mpi_int(1)
    integer, pointer :: v_noext(:) => null()
    integer, pointer :: v_deeq(:) => null()
    integer, pointer :: v_nequ(:) => null()
    integer, pointer :: v_delg(:) => null()
    integer, pointer :: v_owner(:) => null()
    integer, pointer :: v_nulg(:) => null()
    integer, pointer :: v_nuls(:) => null()
    integer, pointer :: v_ddlc(:) => null()
    integer, pointer :: v_nddl(:) => null()
    integer, pointer :: v_gddl(:) => null()
    integer, pointer :: v_tddl(:) => null()
    integer, pointer :: v_nbjo(:) => null()
    integer, pointer :: v_deeg(:) => null()
    integer, pointer :: v_linulg(:) => null()
!
    character(len=8) :: chnbjo
    character(len=8) :: k8bid, mesh, nomgdr
    character(len=19) :: nomlig
    character(len=24) :: owner, nojoie, nojoir, join, linulg
!
!----------------------------------------------------------------------
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
#define zzprno(ili,nunoel,l)  zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
    call jemarq()
    call asmpi_comm('GET', mpicou)
!
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    call jeveuo(numddl//'.NUME.REFN', 'L', jrefn)
    mesh = zk24(jrefn)
    nomgdr = zk24(jrefn + 1)
!
    call jeveuo(mesh//'.DIME', 'L', dime)
    call jeveuo(mesh//'.NULOGL', 'L', vi=v_nulg)
    call jeveuo(mesh//'.NOEX', 'L', vi=v_noext)
!
!   !!! VERIFIER QU'IL N'Y A PAS DE MACRO-ELTS
!   CALCUL DU NOMBRE D'ENTIERS CODES A PARTIR DE LONMAX
    call jelira(jexnum(numddl//'.NUME.PRNO', 1), 'LONMAX', ntot, k8bid)
    nec = ntot/zi(dime) - 2
    call wkvect('&&CRNSTD.NEC', 'V V I', nec, jencod)
    call wkvect('&&CRNSTD.NEC2', 'V V I', nec, jenco2)
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgdr), 'L', iad)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgdr), 'LONMAX', ncmpmx, k8bid)
    call wkvect('&&CRNSTD.CMP', 'V V I', ncmpmx, jcpnec)
    call wkvect('&&CRNSTD.CMP2', 'V V I', ncmpmx, jcpne2)
!
    call jeveuo(numddl//'.NUME.DEEQ', 'L', vi=v_deeq)
    call jeveuo(numddl//'.NUME.NEQU', 'L', vi=v_nequ)
    call jeveuo(numddl//'.NUME.DELG', 'L', vi=v_delg)
    call jeveuo(numddl//'.NUME.NBJO', 'L', vi=v_nbjo)
!
    nbjoin = v_nbjo(1)
    nbddll = v_nequ(1)
!
    call wkvect(numddl//'.NUME.DEEG', 'G V I', 2*nbddll, vi=v_deeg)
    call wkvect(numddl//".NUME.NULS", 'G V I', nbddll, vi=v_nuls)
    v_nuls(:) = -1
!
!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(numddl//'.NUME.PRNO', 'E', idprn1)
    call jeveuo(jexatr(numddl//'.NUME.PRNO', 'LONCUM'), 'L', idprn2)
    call jelira(numddl//'.NUME.PRNO', 'NMAXOC', ntot, k8bid)
!
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbno)
! -- On compte les noeuds locaux proprio
    nbno_lc = 0
    do ino = 1, nbno
        if(v_noext(ino) == rang) then
            nbno_lc = nbno_lc + 1
        end if
    end do
!
! -- Nbr de noeud total
    nbno_gl = nbno_lc
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbno_gl)
!
! -- Nbr de noeud max par proc
    nbno_max = nbno_lc
    call asmpi_comm_vect('MPI_MAX', 'I', sci=nbno_max)
!
! -- On cherche le nombre de composante par noeud physique local proprio
    call wkvect('&&CRNSTD.NODDLL', 'V V I', 2*nbno_max, vi=v_nddl)
    v_nddl(1:2*nbno_max) = -1
    nbno_lc = 0
    numloc = 0
    do ino = 1, nbno
        if(v_noext(ino) == rang) then
            nbno_lc = nbno_lc + 1
! -- Pour les noeuds proprio, on garde le num global et le nombre de ddl du noeud
            v_nddl(2*(nbno_lc-1)+1) = v_nulg(ino)
            v_nddl(2*(nbno_lc-1)+2) = zzprno(1, ino, 2)
            numloc = numloc + zzprno(1, ino, 2)
        end if
    end do
!
    call wkvect('&&CRNSTD.NLGDDL', 'V V I', 2*nbproc*nbno_max, vi=v_gddl)
! -- On envoie tout les ddls - on pourrait optimiser la mémoire
    nbno4 = to_mpi_int(2*nbno_max)
    call asmpi_allgather_i(v_nddl, nbno4, v_gddl, nbno4, mpicou)
!
! -- On trie les ddls par noeuds
    call wkvect('&&CRNSTD.NUMDDL', 'V V I', nbno_gl, vi=v_tddl)
    v_tddl(1:nbno_gl) = -1
    do ino = 1, nbproc*nbno_max
        numero_noeud = v_gddl(2*(ino-1)+1)
        if(numero_noeud .ne. -1) then
! numero_noeud commence à 0
            v_tddl(numero_noeud+1) = v_gddl(2*(ino-1)+2)
        end if
    end do
!
! -- Verif
    do ino = 1, nbno_gl
        ASSERT(v_tddl(ino) >= 0)
    end do
!
! -- On calcule le décalage
    do ino = 2, nbno_gl
        v_tddl(ino) = v_tddl(ino) + v_tddl(ino-1)
    end do
    nbddl_phys_gl = v_tddl(nbno_gl)
    do ino = nbno_gl, 2, -1
        v_tddl(ino) = v_tddl(ino-1)
    end do
    v_tddl(1) = 0
!
! -- On renumerote les noeuds physiques proprio
    call wkvect('&&CRNSTD.DDLLOC', 'V V I', nbno, vi=v_ddlc)
    do i_ddl = 1, nbddll
        numero_noeud = v_deeq((i_ddl-1)*2 + 1)
        numero_cmp   = v_deeq((i_ddl-1)*2 + 2)
        if( numero_noeud.gt.0 .and. numero_cmp.gt.0 ) then
            v_deeg((i_ddl-1)*2 + 1) = v_nulg(numero_noeud) + 1
            v_deeg((i_ddl-1)*2 + 2) = numero_cmp
            if (v_noext(numero_noeud) == rang) then
                v_nuls(i_ddl) = v_tddl(v_nulg(numero_noeud)+1) + v_ddlc(numero_noeud)
                v_ddlc(numero_noeud) = v_ddlc(numero_noeud) + 1
            endif
        endif
    end do
!
! -- On renumerote les noeuds physiques non-proprio
!    Il faut faire comme dans crnlgc.F90
        do i_join = 1, nbjoin
!
!       RECHERCHE DU JOINT
        numpro = v_nbjo(1 + i_join)
        if (numpro .ne. -1) then
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
            num4 = to_mpi_int(i_join)
            numpr4 = to_mpi_int(numpro)
            lgenve1 = nbnoee*(1 + nec) + 1
            lgenvr1 = nbnoer*(1 + nec) + 1
            call wkvect('&&CRNSTD.NOEUD_NEC_E1', 'V V I', lgenvr1, jenvoi1)
            call wkvect('&&CRNSTD.NOEUD_NEC_R1', 'V V I', lgenve1, jrecep1)
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
!
!           On continue si le joint à des DDL
            if(zi(jrecep1) > 0) then
                call wkvect('&&CRNSTD.NUM_DDL_GLOB_E', 'V V I', zi(jrecep1), jenvoi2)
                call wkvect('&&CRNSTD.NUM_DDL_GLOB_R', 'V V I', zi(jenvoi1), jrecep2)
!
                nbddl = 0
                do jaux = 1, nbnoee
                    poscom = (jaux - 1)*(1 + nec) + 1
                    numno1 = zi(jrecep1 + poscom)
!
                    nddl = zzprno(1, numno1, 1)
                    nddlg = v_nuls(nddl)
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
                                nbddl = nbddl + 1
                            endif
                            ico2 = ico2 + 1
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
!
                curpos = 0
                do jaux = 1, nbnoer
                    numno1 = zi(jjoinr + 2*(jaux - 1))
                    nddll = zzprno(1, numno1, 1)
                    nbcmp = zzprno(1, numno1, 2)
                    do icmp = 0, nbcmp - 1
                        ASSERT(zi(jrecep2 + curpos) .ne. -1)
                        v_nuls(nddll + icmp) = zi(jrecep2 + curpos)
                        curpos = curpos + 1
                    enddo
                enddo
                ASSERT(curpos .eq. nb_ddl_envoi)
!
                call jedetr('&&CRNSTD.NUM_DDL_GLOB_E')
                call jedetr('&&CRNSTD.NUM_DDL_GLOB_R')
            endif
!
            call jedetr('&&CRNSTD.NOEUD_NEC_E1')
            call jedetr('&&CRNSTD.NOEUD_NEC_R1')
        endif
    end do
!
! -- On compte les lagranges
    nbddl_lag = 0
    nbddl_lag_gl = 0
    do ili = 2, ntot
        nbno_lili_lc = 0
        call jeexin(jexnum(numddl//'.NUME.PRNO', ili), iret)
        if( iret.ne.0 ) then
            call jelira(jexnum(numddl//'.NUME.PRNO', ili), 'LONMAX', lonmax)
            nbno_prno = lonmax/(nec+2)
            call jenuno(jexnum(numddl//'.NUME.LILI', ili), nomlig)
            owner = nomlig//'.PNOE'
            linulg = nomlig//'.NULG'
            call jeveuo(owner, 'L', vi=v_owner)
            call jeveuo(linulg, 'L', vi=v_linulg)
            do ino = 1, nbno_prno
                ! Le proc est proprio du noeud
                i_ddl = zzprno(ili, ino, 1)
                nbcmp = zzprno(ili, ino, 2)
                ASSERT(nbcmp.eq.1)
                numero_noeud = -nbddl_lag_gl + v_linulg(ino)
                numero_cmp   = v_deeq((i_ddl-1)*2 + 2)
                v_deeg((i_ddl-1)*2 + 1) = numero_noeud
                v_deeg((i_ddl-1)*2 + 2) = numero_cmp
                if( v_owner(ino) == rang ) then
                    nbno_lili_lc = nbno_lili_lc + 1
                    nbddl_lag = nbddl_lag + 1
                    numloc = numloc + 1
                    ASSERT(nbcmp.eq.1)
                    v_nuls(i_ddl) = nbddl_phys_gl - 1 - v_deeg(2*(i_ddl-1)+1)
                endif
            enddo
        endif
!
! -- Nbr de noeud de Lagrange total au ligrel
        nbno_lili_gl = nbno_lili_lc
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nbno_lili_gl)
        nbddl_lag_gl = nbddl_lag_gl + nbno_lili_gl
    enddo
!
!   Nombre total de Lagrange
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbddl_lag)
    ASSERT(nbddl_lag == nbddl_lag_gl)
!
!   Nombre total de degré de liberté
    nbddll_gl = numloc
    call asmpi_comm_vect('MPI_SUM', 'I', sci=nbddll_gl)
!
! -- Verif nb ddl total
    ASSERT(nbddll_gl == nbddl_phys_gl + nbddl_lag_gl)
!
! -- On complete avec les joins
    do ili = 2, ntot
        call jeexin(jexnum(numddl//'.NUME.PRNO', ili), iret)
        if( iret.ne.0 ) then
            call jenuno(jexnum(numddl//'.NUME.LILI', ili), nomlig)
            join = nomlig//".NBJO"
            call jeexin(join, iret)
            if( iret.eq.0 ) cycle
            call jeveuo(join, 'L', jjoin)
            call jelira(join, 'LONMAX', nbjoin)
            do i_join = 1, nbjoin
                numpro = zi(jjoin+i_join-1)
                if( numpro.ne.-1 ) then
                    numpr4 = to_mpi_int(numpro)
                    num4 = to_mpi_int(i_join)
                    call codent(numpro, 'G', chnbjo)
                    nojoie = nomlig//'.E'//chnbjo
                    nojoir = nomlig//'.R'//chnbjo

                    call jeexin(nojoie, iret1)
                    if( iret1.ne.0 ) then
                        call jeveuo(nojoie, 'L', jjoine)
                        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
                        call wkvect('&&CRNSTD.NUM_DDL_GLOB_E', 'V V I', nbnoee, jnujoi1)
                        do jaux = 1, nbnoee
                            numnoe = -zi(jjoine+jaux-1)
                            ASSERT(zzprno(ili, numnoe, 2).eq.1)
                            nddll = zzprno(ili, numnoe, 1)
                            zi(jnujoi1+jaux-1) = v_nuls(nddll)
                        enddo
                        n4e = nbnoee
                    endif

                    call jeexin(nojoir, iret2)
                    if( iret2.ne.0 ) then
                        call jeveuo(nojoir, 'L', jjoinr)
                        call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
                        call wkvect('&&CRNSTD.NUM_DDL_GLOB_R', 'V V I', nbnoer, jnujoi2)
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
                            v_nuls(nddll) = zi(jnujoi2+jaux-1)
                        enddo
                    endif
                    call jedetr('&&CRNSTD.NUM_DDL_GLOB_E')
                    call jedetr('&&CRNSTD.NUM_DDL_GLOB_R')
                endif
            enddo
        endif
    enddo
!
! -- Verif finale
    do i_ddl = 1, nbddll
        ASSERT(v_nuls(i_ddl) .ne. -1)
    end do
!
! -- On détruit
    call jedetr('&&CRNSTD.NLGDDL')
    call jedetr('&&CRNSTD.NODDLL')
    call jedetr('&&CRNSTD.NUMDLL')
    call jedetr('&&CRNSTD.DDLLOC')
    call jedetr('&&CRNSTD.NEC')
    call jedetr('&&CRNSTD.NEC2')
    call jedetr('&&CRNSTD.CMP')
    call jedetr('&&CRNSTD.CMP2')
!
    call jedema()
#else
    character(len=14) :: k14
    k14 = numddl
#endif
!
end subroutine
