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

subroutine macr78(nomres, trange, typres)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdgeph.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rstran.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: nomres
    character(len=16) :: typres
    character(len=19) :: trange
! IN  : NOMRES : NOM UTILISATEUR POUR LA COMMANDE REST_COND_TRAN
! IN  : TYPRES : TYPE DE RESULTAT : 'DYNA_TRANS'
! IN  : TRANGE : NOM UTILISATEUR DU CONCEPT TRAN_GENE AMONT
!       NOMCMD : NOM DE LA COMMANDE : 'REST_COND_TRAN'
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    complex(kind=8) :: cbid
    character(len=8) :: k8b, basemo, mailla, nomin, nomcmp(6), macrel, lintf
    character(len=8) :: nomnol, nogdsi, maya
    character(len=14) :: numddl
    character(len=16) :: concep, champ(8)
    character(len=19) :: kinst, knume, cham19
    character(len=24) :: chamno, nomcha, numedd, nprno
!      CHARACTER*3  TREDU
    aster_logical :: lredu
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iaprno, iarc0, iarch
    integer(kind=8) :: ibid, icmp, iddl, im
    integer(kind=8) :: inoe, inu0, inum, iret, iretou, ivale, j
    integer(kind=8) :: jinst, jnume, k, ldnew
    integer(kind=8) :: linst, lnocm2, lnocmp, lpa2, lpar, n0
    integer(kind=8) :: n1, nbcham, nbec, nbinst, nbmdef, nbmdyn, nbmode
    integer(kind=8) :: nbndef, nbndyn, nbnoe, nbntot, nbtdyn, nec, neq
    integer(kind=8) :: nes, nmc, tmod(1)
    real(kind=8) :: rbid
    real(kind=8), pointer :: base(:) => null()
    character(len=8), pointer :: noecmp(:) => null()
    real(kind=8), pointer :: restr(:) => null()
    integer(kind=8), pointer :: lino(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
    character(len=24), pointer :: mael_refe(:) => null()
!-----------------------------------------------------------------------
    data nomcmp/'DX      ', 'DY      ', 'DZ      ',&
     &               'DRX     ', 'DRY     ', 'DRZ     '/
!     ------------------------------------------------------------------
    call jemarq()
    nomin = trange(1:8)
!      TYPRES = 'DYNA_TRANS'
    call gettco(nomin, concep)
!
!     --- RECUPERATION DES ENTITES DU MAILLAGE SUR LESQUELLES ---
!     ---                PORTE LA RESTITUTION                 ---
!      TOUSNO = .TRUE.
!
    lredu = .false.
!      CALL GETVTX ( ' ', 'REDUC' , 1,IARG,1, TREDU, N2 )
!      IF (TREDU.EQ.'OUI') LREDU = .TRUE.
    call getvid(' ', 'MACR_ELEM_DYNA', scal=macrel, nbret=nmc)
    call jeveuo(macrel//'.MAEL_REFE', 'L', vk24=mael_refe)
    basemo = mael_refe(1) (1:8)
    call rsorac(basemo, 'LONUTI', 0, rbid, k8b, &
                cbid, rbid, k8b, tmod, 1, &
                ibid)
    nbmode = tmod(1)
    call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=numedd)
    call dismoi('NOM_MAILLA', numedd(1:14), 'NUME_DDL', repk=mailla)
    call dismoi('REF_INTD_PREM', basemo, 'RESU_DYNA', repk=lintf)
    call jelira(jexnum(lintf//'.IDC_LINO', 1), 'LONMAX', nbnoe)
    call dismoi('NB_MODES_STA', basemo, 'RESULTAT', repi=nbmdef)
!      CALL BMNBMD(BASEMO,'DEFORMEE',NBMDEF)
    nbmdyn = nbmode-nbmdef
    call jelira(macrel//'.LINO', 'LONMAX', nbntot)
    nec = nbmode/nbntot
    nbndyn = nbmdyn/nec
    nbndef = nbntot-nbndyn
!      NBNDE2 = NBMDEF/NEC
!      ASSERT(NBNDEF.EQ.NBNDE2)
    if (nbmdef .ne. 0) then
        call rsadpa(basemo, 'L', 1, 'NOEUD_CMP', nbmdyn+1, &
                    0, sjv=lnocmp, styp=k8b)
        if (zk16(lnocmp) .eq. ' ') then
            lredu = .true.
            nec = nbmode/nbntot
            nbndyn = nbmdyn/nec
            nbndef = nbntot-nbndyn
        else
            k = 1
31          continue
            if ((k+1) .gt. nbmdef) then
                nes = k
                goto 32
            end if
            call rsadpa(basemo, 'L', 1, 'NOEUD_CMP', nbmdyn+k+1, &
                        0, sjv=lnocm2, styp=k8b)
            if (zk16(lnocmp) (1:8) .ne. zk16(lnocm2) (1:8)) then
                nes = k
                goto 32
            else
                k = k+1
                goto 31
            end if
32          continue
            nbndef = nbmdef/nes
            nbndyn = nbntot-nbndef
            if (nbmdyn .ne. 0) then
                nec = nbmdyn/nbndyn
            else
                nec = nes
            end if
        end if
    end if
!       CREATION DU TABLEAU NOEUD-COMPOSANTE ASSOCIES AUX MODES
    AS_ALLOCATE(vk8=noecmp, size=2*nbmode)
    call jeveuo(macrel//'.LINO', 'L', vi=lino)
    if (lredu) then
        nbtdyn = nbntot
    else
        nbtdyn = nbndyn
        do i = nbmdyn+1, nbmode
            call rsadpa(basemo, 'L', 1, 'NOEUD_CMP', i, &
                        0, sjv=lnocmp, styp=k8b)
            noecmp(1+2*i-2) = zk16(lnocmp) (1:8)
            noecmp(1+2*i-1) = zk16(lnocmp) (9:16)
        end do
    end if
!
    do i = 1, nbtdyn
        nomnol = int_to_char8(lino(i))
        do j = 1, nec
            noecmp(1+2*nec*(i-1)+2*j-2) = nomnol
            noecmp(1+2*nec*(i-1)+2*j-1) = nomcmp(j)
        end do
    end do
!        CALL GETVID(' ','NUME_DDL',1,IARG,1,K8B,IBID)
!        IF (IBID.NE.0) THEN
!          CALL GETVID(' ','NUME_DDL',1,1,1,NUMEDD,IBID)
!          NUMEDD = NUMEDD(1:14)//'.NUME'
!        ENDIF
    numddl = numedd(1:14)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
    AS_ALLOCATE(vr=base, size=nbmode*neq)
    call copmod(basemo, bmodr=base, numer=numddl)
    call getvtx(' ', 'TOUT_CHAM', nbval=0, nbret=n0)
    if (n0 .ne. 0) then
        nbcham = 3
        champ(1) = 'DEPL'
        champ(2) = 'VITE'
        champ(3) = 'ACCE'
    else
        call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbcham = -n1
            call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=champ, nbret=n1)
        else
            call utmess('A', 'ALGORITH10_93')
            goto 999
        end if
    end if
    knume = '&&MACR78.NUM_RANG'
    kinst = '&&MACR78.INSTANT'
    call rstran('NON', trange, ' ', 1, kinst, &
                knume, nbinst, iretou)
    if (iretou .ne. 0) then
        call utmess('F', 'UTILITAI4_24')
    end if
    call jeexin(kinst, iret)
    if (iret .gt. 0) then
        call jeveuo(kinst, 'L', jinst)
        call jeveuo(knume, 'L', jnume)
    end if
!
    call jeexin(trange//'.ORDR', iret)
    if (iret .ne. 0) then
        call jeveuo(trange//'.ORDR', 'L', vi=ordr)
!
    end if
!
!
!     --- CREATION DE LA SD RESULTAT ---
    call rscrsd('G', nomres, typres, nbinst)
!
    AS_ALLOCATE(vr=restr, size=nbmode)
    call rsexch('F', nomin, 'DEPL', 1, cham19, &
                iret)
    call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=maya)
! ATTENTION MAILLA MAILLAGE DU MACRO_ELEM DE RESTITUTION
!  ET MAYA MAILLAGE DU RESULTAT SUR MODELE SOUS-STRUC-STAT
    call dismoi('NOM_GD', cham19, 'CHAMP', repk=nogdsi)
    call dismoi('NB_EC', nogdsi, 'GRANDEUR', repi=nbec)
!
    call dismoi('NUME_EQUA', cham19, 'CHAMP', repk=nprno)
    nprno = nprno(1:19)//'.PRNO'
    call jeveuo(jexnum(nprno, 1), 'L', iaprno)
    do i = 1, nbcham
        do iarc0 = 1, nbinst
            inu0 = zi(jnume+iarc0-1)
            inum = ordr(inu0)
            iarch = iarc0-1
            call rsexch('F', nomin, champ(i) (1:4), inum, nomcha, &
                        iret)
            nomcha = nomcha(1:19)//'.VALE'
            call jeveuo(nomcha, 'L', ivale)
            do im = 1, nbmode
                nomnol = noecmp(1+2*im-2)
                inoe = char8_to_int(nomnol)
                if (noecmp(1+2*im-1) .eq. 'DX') icmp = 1
                if (noecmp(1+2*im-1) .eq. 'DY') icmp = 2
                if (noecmp(1+2*im-1) .eq. 'DZ') icmp = 3
                if (noecmp(1+2*im-1) .eq. 'DRX') icmp = 4
                if (noecmp(1+2*im-1) .eq. 'DRY') icmp = 5
                if (noecmp(1+2*im-1) .eq. 'DRZ') icmp = 6
                iddl = zi(iaprno-1+(nbec+2)*(inoe-1)+1)
                restr(im) = zr(ivale+iddl-1+icmp-1)
            end do
            call rsexch(' ', nomres, champ(i) (1:4), iarch, chamno, &
                        iret)
            call vtcreb(chamno, 'G', 'R', &
                        nume_ddlz=numedd, &
                        nb_equa_outz=neq)
            call jeveuo(chamno(1:19)//'.VALE', 'E', ldnew)
            call mdgeph(neq, nbmode, base, restr, zr(ldnew))
            call rsnoch(nomres, champ(i) (1:4), iarch)
            if (i .eq. 1) then
                call rsadpa(nomres, 'E', 1, 'INST', iarch, &
                            0, sjv=linst, styp=k8b)
                zr(linst) = zr(jinst+iarc0-1)
                call rsadpa(nomin, 'L', 1, 'MODELE', inum, &
                            0, sjv=lpa2, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'MODELE', iarch, &
                            0, sjv=lpar, styp=k8b)
                zk8(lpar) = zk8(lpa2)
                call rsadpa(nomin, 'L', 1, 'CHAMPMAT', inum, &
                            0, sjv=lpa2, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'CHAMPMAT', iarch, &
                            0, sjv=lpar, styp=k8b)
                zk8(lpar) = zk8(lpa2)
                call rsadpa(nomin, 'L', 1, 'CARAELEM', inum, &
                            0, sjv=lpa2, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'CARAELEM', iarch, &
                            0, sjv=lpar, styp=k8b)
                zk8(lpar) = zk8(lpa2)
            end if
        end do
    end do
!
    call refdcp(basemo, nomres)
!
!
! --- MENAGE
!
    AS_DEALLOCATE(vk8=noecmp)
    AS_DEALLOCATE(vr=base)
    call jedetr('&&MACR78.NUM_RANG')
    call jedetr('&&MACR78.INSTANT')
    AS_DEALLOCATE(vr=restr)
!
    call titre()
999 continue
!
    call jedema()
end subroutine
