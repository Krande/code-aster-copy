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

subroutine rfhge2(harmge)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/mdgep5.h"
#include "asterfort/posddl.h"
#include "asterfort/rstran.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zxtrac.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: harmge
!
!     OPERATEUR "RECU_FONCTION"  MOT CLE "HARM_GENE"
!     ------------------------------------------------------------------
    character(len=4) :: interp(2)
    character(len=8) :: crit, noeud, cmp, noma, basemo
    character(len=8) :: intres
    character(len=14) :: nume
    character(len=16) :: nomcmd, typcon, nomcha
    character(len=19) :: nomfon, knume, kinst, resu
    character(len=24) :: nogno, valk(2)
    complex(kind=8) :: crep, cbid
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iagno, idbase, iddl
    integer(kind=8) :: ie, ierd, ign2, ii, ino, inoeud, iordr
    integer(kind=8) :: iret, itresu, jinst, jj, lfon, lg1, lg2
    integer(kind=8) :: lordr, lpro, lvar, n1, n2
    integer(kind=8) :: n3, nbinsg, nbmode, nbordr
    integer(kind=8) :: neq, ngn, numcmp
    real(kind=8) :: epsi
    complex(kind=8), pointer :: vectgene(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    call jemarq()
!
    call getres(nomfon, typcon, nomcmd)
!
    resu = harmge
    interp(1) = 'LIN '
    interp(2) = 'LIN '
    intres = 'NON '
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERP_NUME', scal=intres, nbret=n1)
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n1)
    if (n1 .eq. 1) interp(2) = interp(1)
!
    noeud = ' '
    cmp = ' '
    call getvtx(' ', 'NOM_CMP', scal=cmp, nbret=n1)
    call getvtx(' ', 'NOM_CHAM', scal=nomcha, nbret=n2)
    call getvtx(' ', 'NOEUD', scal=noeud, nbret=n3)
    call getvtx(' ', 'GROUP_NO', scal=nogno, nbret=ngn)
!
    call jeexin(resu//'.'//nomcha(1:4), iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_23', sk=nomcha)
    end if
    call jeveuo(resu//'.'//nomcha(1:4), 'L', itresu)
!
    knume = '&&RFHGE2.NUME_ORDR'
    kinst = '&&RFHGE2.FREQUENCE'
    call rstran(intres, resu, ' ', 1, kinst, &
                knume, nbordr, ie)
    if (ie .ne. 0) then
        call utmess('F', 'UTILITAI4_15')
    end if
    call jeexin(kinst, iret)
    if (iret .gt. 0) then
        call jeveuo(kinst, 'L', jinst)
        call jeveuo(knume, 'L', lordr)
    end if
!
!     --- CREATION DE LA FONCTION ---
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
    zk24(lpro) = 'FONCT_C         '
    zk24(lpro+1) = interp(1)//interp(2)
    zk24(lpro+2) = 'FREQ            '
    zk24(lpro+3) = nomcha
    zk24(lpro+4) = 'EE              '
    zk24(lpro+5) = nomfon(1:19)
!
! --- LA FONCTION EST LA CONCATENATION DE DEUX VECTEURS:
! --- ABSCISSES +  ( PARTIE REELLE | PARTIE IMAGINAIRE )
    call wkvect(nomfon//'.VALE', 'G V R', 3*nbordr, lvar)
!
    call jeveuo(resu//'.DESC', 'L', vi=desc)
    nbmode = desc(2)
    call getvis(' ', 'NUME_CMP_GENE', scal=numcmp, nbret=n1)
    lfon = lvar+nbordr
!
! --- CAS OU D'UNE VARIABLE GENERALISEE
!
    if (n1 .ne. 0) then
        if (numcmp .gt. nbmode) then
            call utmess('F', 'UTILITAI4_14')
        end if
!
        jj = 0
        if (intres(1:3) .ne. 'NON') then
! ---   CAS OU ON INTERPOLE
            call utmess('E', 'ALGORITH11_79')
        else
! ---   CAS OU ON N'INTERPOLE PAS
            do iordr = 0, nbordr-1
                ii = zi(lordr+iordr)
                zr(lvar+iordr) = zr(jinst+iordr)
                crep = zc(itresu+nbmode*(ii-1)+numcmp-1)
                zr(lfon+jj) = dble(crep)
                jj = jj+1
                zr(lfon+jj) = dimag(crep)
                jj = jj+1
            end do
        end if
    else
!
! --- CAS D'UNE VARIABLE PHYSIQUE
!
        call dismoi('BASE_MODALE', resu, 'RESU_DYNA', repk=basemo)
        call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=nume)
        call dismoi('NOM_MAILLA', nume, 'NUME_DDL', repk=noma)
!
! ---   RECUPERATION DE LA BASE MODALE DANS UN VECTEUR DE TRAVAIL
        call dismoi('NB_EQUA', nume, 'NUME_DDL', repi=neq)
        call wkvect('&&RFHGE2.VECT.PROPRE', 'V V R', neq*nbmode, idbase)
        call copmod(basemo, numer=nume, bmodr=zr(idbase))
!
! --- TRAITEMENT D'UN GROUP DE NOEUDS SEUELEMENT
        if (ngn .ne. 0) then
            call jenonu(jexnom(noma//'.GROUPENO', nogno), ign2)
            if (ign2 .le. 0) then
                call utmess('F', 'ELEMENTS_67', sk=nogno)
            end if
            call jeveuo(jexnum(noma//'.GROUPENO', ign2), 'L', iagno)
            ino = zi(iagno)
            noeud = int_to_char8(ino)
        end if
        call posddl('NUME_DDL', nume, noeud, cmp, inoeud, &
                    iddl)
        if (inoeud .eq. 0) then
            lg1 = lxlgut(noeud)
            call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
        else if (iddl .eq. 0) then
            lg1 = lxlgut(noeud)
            lg2 = lxlgut(cmp)
            valk(1) = cmp(1:lg2)
            valk(2) = noeud(1:lg1)
            call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
        end if
! --- INTERPOLATION PROPREMENT DITE (ESPACE PHYSIQUE)
        jj = 0
        if (intres(1:3) .ne. 'NON') then
! ---   CAS OU ON INTERPOLE
            call jeveuo(resu//'.DISC', 'L', vr=disc)
            call jelira(resu//'.DISC', 'LONMAX', nbinsg)
            AS_ALLOCATE(vc=vectgene, size=nbmode)
            do iordr = 0, nbordr-1
!             EXTRACTION ET INTERPOLATION
                call zxtrac(intres, epsi, crit, nbinsg, disc, &
                            zr(jinst+iordr), zc(itresu), nbmode, vectgene, ierd)
!             PASSAGE EN BASE PHYSIQUE
                call mdgep5(neq, nbmode, zr(idbase), vectgene, iddl, &
                            crep)
!             REMPLISSAGE DES TROIS VECTEURS DE LA FONCTION
                zr(lvar+iordr) = zr(jinst+iordr)
                zr(lfon+jj) = dble(crep)
                jj = jj+1
                zr(lfon+jj) = dimag(crep)
                jj = jj+1
            end do
            AS_DEALLOCATE(vc=vectgene)
!
        else
! ---   CAS OU ON N'INTERPOLE PAS
            do iordr = 0, nbordr-1
                ii = zi(lordr+iordr)
!             PASSAGE EN BASE PHYSIQUE
                call mdgep5(neq, nbmode, zr(idbase), zc(itresu+nbmode*(ii-1)), iddl, &
                            crep)
                zr(lvar+iordr) = zr(jinst+iordr)
                zr(lfon+jj) = dble(crep)
                jj = jj+1
                zr(lfon+jj) = dimag(crep)
                jj = jj+1
            end do
!
        end if
    end if
!
    call jedetr('&&RFHGE2.VECT.PROPRE')
!
!     ---------------------------------------------------------------
    call jedetr(knume)
    call jedetr(kinst)
!
    call jedema()
end subroutine
