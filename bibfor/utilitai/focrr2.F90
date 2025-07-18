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

subroutine focrr2(nomfon, resu, base, nomcha, maille, &
                  noeud, cmp, npoint, nusp, ivari, nomvari, &
                  ier)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/focrrs.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/posddl.h"
#include "asterfort/rsbary.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslipa.h"
#include "asterfort/rsutro.h"
#include "asterfort/utch19.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=1), intent(in) :: base
    character(len=8), intent(in) :: maille, noeud, cmp
    character(len=16), intent(in) :: nomcha, nomvari
    character(len=19), intent(in) :: nomfon, resu
!     RECUPERATION D'UNE FONCTION DANS UNE STRUCTURE "RESULTAT"
!     ------------------------------------------------------------------
! VAR : NOMFON : NOM DE LA FONCTION
! IN  : RESU   : NOM DE LA STRUCTURE RESULTAT
! IN  : BASE   : BASE OU L'ON CREE LA FONCTION
! IN  : NOMCHA : NOM DU CHAMP
! IN  : NOEUD  : NOEUD
! IN  : MAILLE : MAILE
! IN  : CMP    : COMPOSANTE
! IN  : NPOINT : NUMERO DU POINT ( CAS DES CHAM_ELEMS )
! IN  : NUSP   : NUMERO DU SOUS-POINT ( CAS DES CHAM_ELEMS )
! IN  : IVARI   : NUMERO DE LA CMP (POUR VARI_R)
! IN  : NOMVARI : NOM DE LA CMP (POUR VARI_R)
! OUT : IER    : CODE RETOUR, = 0 : OK
!     ------------------------------------------------------------------
    character(len=1) :: type
    character(len=4) :: typch2
    character(len=8) :: k8b, noma, nogd, nomacc, noeuz
    character(len=16) :: nomcmd, typcon, typcha, typres
    character(len=19) :: listr, profch, profc2, ch1, ch2
    character(len=24) :: valk(2)
    real(kind=8) :: dimag
    real(kind=8) :: valr(2)
    complex(kind=8) :: valc1, valc2
    integer(kind=8) :: npoinz, nuspz
    integer(kind=8) :: i, i1, i2, iacces, iatach, iatava
    integer(kind=8) :: ibid, iddl1, iddl2, ier, ierd, ierr1
    integer(kind=8) :: ierr2, ii, inoeud, iordr, ip1, ip2, iposit
    integer(kind=8) :: iret, ivari, jinst, jlir8, l1, l2, lfon
    integer(kind=8) :: lg1, lg2, lpro, lval1, lvar
    integer(kind=8) :: n1, n2, n3, n4, nbinst, nbordr, npoint
    integer(kind=8) :: nusp, vali1, vali2
    real(kind=8) :: r1, r2, rbase, rval, valr1, valr2
    aster_logical, pointer :: lexi(:) => null()
    real(kind=8), pointer :: nlval2(:) => null()
    real(kind=8), parameter :: prec = 1.0d-10
    character(len=8), parameter :: crit = 'ABSOLU'
!
!-----------------------------------------------------------------------
    call jemarq()
!
    ier = 0
    call getres(k8b, typcon, nomcmd)
    call gettco(resu, typres)
!
    call getvr8(' ', 'INST', nbval=0, nbret=n1)
    call getvid(' ', 'LIST_INST', nbval=0, nbret=n2)
    call getvr8(' ', 'FREQ', nbval=0, nbret=n3)
    call getvid(' ', 'LIST_FREQ', nbval=0, nbret=n4)
!
    if (typres(1:10) .eq. 'DYNA_HARMO') then
        nomacc = 'FREQ'
        if (n1+n2 .ne. 0) then
            call utmess('F', 'UTILITAI_95')
        end if
        if (n3+n4 .eq. 0) then
            call focrrs(nomfon, resu, base, nomcha, maille, &
                        noeud, cmp, npoint, nusp, ivari, nomvari, &
                        ier)
            goto 40
        end if
        if (n3 .ne. 0) then
            nbinst = -n3
            call wkvect('&&FOCRR2.INST', 'V V R', nbinst, jinst)
            call getvr8(' ', 'FREQ', nbval=nbinst, vect=zr(jinst), nbret=n3)
        else
            call getvid(' ', 'LIST_FREQ', scal=listr, nbret=n4)
            call jeveuo(listr//'.VALE', 'L', jinst)
            call jelira(listr//'.VALE', 'LONMAX', nbinst)
        end if
    else
        nomacc = 'INST'
        if (n3+n4 .ne. 0) then
            call utmess('F', 'UTILITAI_96')
        end if
        if (n1+n2 .eq. 0) then
            call focrrs(nomfon, resu, base, nomcha, maille, &
                        noeud, cmp, npoint, nusp, ivari, nomvari, &
                        ier)
            goto 40
        end if
        if (n1 .ne. 0) then
            nbinst = -n1
            call wkvect('&&FOCRR2.INST', 'V V R', nbinst, jinst)
            call getvr8(' ', 'INST', nbval=nbinst, vect=zr(jinst), nbret=n1)
        else
            call getvid(' ', 'LIST_INST', scal=listr, nbret=n2)
            call jeveuo(listr//'.VALE', 'L', jinst)
            call jelira(listr//'.VALE', 'LONMAX', nbinst)
        end if
    end if
!
!     --- REMPLISSAGE DU .PROL ---
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', base//' V K24', 6, lpro)
    if (typres(1:10) .eq. 'DYNA_HARMO') then
        zk24(lpro) = 'FONCT_C'
    else
        zk24(lpro) = 'FONCTION'
    end if
    zk24(lpro+1) = 'LIN LIN '
    zk24(lpro+2) = nomacc
    zk24(lpro+3) = cmp
    zk24(lpro+4) = 'EE      '
    zk24(lpro+5) = nomfon
!
    if (typres(1:10) .eq. 'DYNA_HARMO') then
        call wkvect(nomfon//'.VALE', base//' V R', 3*nbinst, lvar)
    else
        call wkvect(nomfon//'.VALE', base//' V R', 2*nbinst, lvar)
    end if
    lfon = lvar+nbinst
!
    call jenonu(jexnom(resu//'.NOVA', nomacc), iacces)
    ASSERT(iacces .ne. 0)
    call jeveuo(jexnum(resu//'.TAVA', iacces), 'L', iatava)
!
    call rslipa(resu, nomacc, '&&FOCRR2.LIR8', jlir8, nbordr)
!
!
!     -- ON REPERE QUELS SONT LES CHAMPS EXISTANT REELLEMENT:
    AS_ALLOCATE(vl=lexi, size=nbordr)
    call jenonu(jexnom(resu//'.DESC', nomcha), ibid)
    call jeveuo(jexnum(resu//'.TACH', ibid), 'L', iatach)
    do i = 1, nbordr
        if (zk24(iatach-1+i) (1:1) .eq. ' ') then
            lexi(i) = .false.
        else
            lexi(i) = .true.
        end if
    end do
!
    rval = zr(jinst)
    call rsbary(zr(jlir8), nbordr, ASTER_FALSE, lexi, rval, &
                i1, i2, iposit, prec, crit)
    call rsutro(resu, i1, ip1, ierr1)
    call rsexch('F', resu, nomcha, ip1, ch1, &
                ierd)
    call dismoi('TYPE_SUPERVIS', ch1, 'CHAMP', repk=typcha)
!
!               ----- EXTRACTION SUR UN "CHAM_NO" -----
!
    if (typcha(1:7) .eq. 'CHAM_NO') then
        call dismoi('NUME_EQUA', ch1, 'CHAM_NO', repk=profch)
        call dismoi('NOM_MAILLA', ch1, 'CHAM_NO', repk=noma)
        call posddl('CHAM_NO', ch1, noeud, cmp, inoeud, &
                    iddl1)
        if (inoeud .eq. 0) then
            lg1 = lxlgut(noeud)
            call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
        else if (iddl1 .eq. 0) then
            lg1 = lxlgut(noeud)
            lg2 = lxlgut(cmp)
            valk(1) = cmp(1:lg2)
            valk(2) = noeud(1:lg1)
            call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
        end if
        iddl2 = iddl1
        do iordr = 0, nbinst-1
            call jemarq()
!
            rval = zr(jinst+iordr)
            call rsbary(zr(jlir8), nbordr, ASTER_FALSE, lexi, rval, &
                        i1, i2, iposit, prec, crit)
            if (iposit .eq. -2) then
                valr(1) = rval
                call utmess('F', 'UTILITAI6_16', sr=valr(1))
!
!           -- PROLONGEMENT A GAUCHE:
!           -------------------------
            else if (iposit .eq. -1) then
                valr(1) = rval
                valr(2) = zr(jlir8)
                call utmess('F', 'UTILITAI6_17', nr=2, valr=valr)
!
!           -- PROLONGEMENT A DROITE:
!           -------------------------
            else if (iposit .eq. 1) then
                valr(1) = rval
                valr(2) = zr(jlir8+nbordr-1)
                call utmess('F', 'UTILITAI6_18', nr=2, valr=valr)
            end if
!
            call rsutro(resu, i1, ip1, ierr1)
            call rsutro(resu, i2, ip2, ierr2)
            ASSERT(ierr1+ierr2 .le. 0)
            rbase = zr(jlir8-1+i2)-zr(jlir8-1+i1)
!
            call rsexch('F', resu, nomcha, ip1, ch1, &
                        l1)
            call rsexch('F', resu, nomcha, ip2, ch2, &
                        l2)
!
            call dismoi('NUME_EQUA', ch1, 'CHAM_NO', repk=profc2)
            if (profc2 .ne. profch) then
                profch = profc2
                call posddl('CHAM_NO', ch1, noeud, cmp, inoeud, &
                            iddl1)
                if (inoeud .eq. 0) then
                    lg1 = lxlgut(noeud)
                    call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
                else if (iddl1 .eq. 0) then
                    lg1 = lxlgut(noeud)
                    lg2 = lxlgut(cmp)
                    valk(1) = cmp(1:lg2)
                    valk(2) = noeud(1:lg1)
                    call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
                end if
                iddl2 = iddl1
            end if
!
            if (rbase .eq. 0.0d0) then
                call jeveuo(ch1//'.VALE', 'L', lval1)
                zr(lvar+iordr) = rval
                zr(lfon+iordr) = zr(lval1+iddl1-1)
                goto 22
            end if
            r1 = (zr(jlir8-1+i2)-rval)/rbase
            r2 = (rval-zr(jlir8-1+i1))/rbase
!
            call dismoi('NUME_EQUA', ch2, 'CHAM_NO', repk=profc2)
            if (profc2 .ne. profch) then
                profch = profc2
                call posddl('CHAM_NO', ch2, noeud, cmp, inoeud, &
                            iddl2)
                if (inoeud .eq. 0) then
                    lg1 = lxlgut(noeud)
                    call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
                else if (iddl2 .eq. 0) then
                    lg1 = lxlgut(noeud)
                    lg2 = lxlgut(cmp)
                    valk(1) = cmp(1:lg2)
                    valk(2) = noeud(1:lg1)
                    call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
                end if
            end if
!
            call jeveuo(ch1//'.VALE', 'L', lval1)
            call jeveuo(ch2//'.VALE', 'L', vr=nlval2)
            zr(lvar+iordr) = rval
            zr(lfon+iordr) = r1*zr(lval1+iddl1-1)+r2*nlval2(iddl2)
!
            iddl1 = iddl2
22          continue
            call jedema()
        end do
!
!               ----- EXTRACTION SUR UN "CHAM_ELEM" -----
!
    else if (typcha(1:9) .eq. 'CHAM_ELEM') then
        noeuz = noeud
        npoinz = npoint
        nuspz = nusp
! ---    VERIFICATION DE LA PRESENCE DES MOTS CLE GROUP_MA (OU MAILLE)
! ---    ET GROUP_NO (OU NOEUD OU POINT) DANS LE CAS D'UN CHAM_ELEM
!        -------------------------------------------------------------
        call dismoi('TYPE_CHAMP', ch1, 'CHAMP', repk=typch2)
        if (typch2 .eq. 'ELEM') then
            npoinz = 1
            nuspz = 1
            noeuz = ' '
            if (maille .eq. ' ') then
                call utmess('F', 'CHAMPS_11')
            end if
        else if (typch2 .eq. 'ELNO') then
            nuspz = 1
            if (maille .eq. ' ' .or. (noeud .eq. ' ' .and. npoint .eq. 0)) then
                call utmess('F', 'CHAMPS_12')
            end if
        else
            if (maille .eq. ' ' .or. npoint .eq. 0) then
                call utmess('F', 'CHAMPS_13')
            end if
        end if
        call dismoi('NOM_MAILLA', ch1, 'CHAM_ELEM', repk=noma)
        call dismoi('NOM_GD', ch1, 'CHAM_ELEM', repk=nogd)
        call dismoi('TYPE_SCA', nogd, 'GRANDEUR', repk=type)
!
        ii = 0
        do iordr = 0, nbinst-1
            call jemarq()
!
            rval = zr(jinst+iordr)
            call rsbary(zr(jlir8), nbordr, ASTER_FALSE, lexi, rval, &
                        i1, i2, iposit, prec, crit)
            if (iposit .eq. -2) then
                valr(1) = rval
                call utmess('F', 'UTILITAI6_16', sr=valr(1))
!
!           -- PROLONGEMENT A GAUCHE:
!           -------------------------
            else if (iposit .eq. -1) then
                valr(1) = rval
                valr(2) = zr(jlir8)
                call utmess('F', 'UTILITAI6_17', nr=2, valr=valr)
!
!           -- PROLONGEMENT A DROITE:
!           -------------------------
            else if (iposit .eq. 1) then
                valr(1) = rval
                valr(2) = zr(jlir8+nbordr-1)
                call utmess('F', 'UTILITAI6_18', nr=2, valr=valr)
            end if
!
            call rsutro(resu, i1, ip1, ierr1)
            call rsutro(resu, i2, ip2, ierr2)
            ASSERT(ierr1+ierr2 .le. 0)
            rbase = zr(jlir8-1+i2)-zr(jlir8-1+i1)
!
            call rsexch('F', resu, nomcha, ip1, ch1, &
                        l1)
            call rsexch('F', resu, nomcha, ip2, ch2, &
                        l2)
!
            if (rbase .eq. 0.0d0) then
                call utch19(ch1, noma, maille, noeuz, npoinz, &
                            nuspz, ivari, cmp, type, valr1, &
                            valc1, vali1, iret)
                ASSERT(iret .eq. 0)
                zr(lvar+iordr) = rval
                if (type .eq. 'R') then
                    zr(lfon+iordr) = valr1
                else
                    zr(lfon+ii) = dble(valc1)
                    ii = ii+1
                    zr(lfon+ii) = dimag(valc1)
                    ii = ii+1
                end if
                goto 32
            end if
            r1 = (zr(jlir8-1+i2)-rval)/rbase
            r2 = (rval-zr(jlir8-1+i1))/rbase
!
            call utch19(ch1, noma, maille, noeuz, npoinz, &
                        nuspz, ivari, cmp, type, valr1, &
                        valc1, vali1, iret)
            ASSERT(iret .eq. 0)
            call utch19(ch2, noma, maille, noeuz, npoinz, &
                        nuspz, ivari, cmp, type, valr2, &
                        valc2, vali2, iret)
            ASSERT(iret .eq. 0)
!
            zr(lvar+iordr) = rval
            if (type .eq. 'R') then
                zr(lfon+iordr) = r1*valr1+r2*valr2
            else
                zr(lfon+ii) = dble(r1*valc1+r2*valc2)
                ii = ii+1
                zr(lfon+ii) = dimag(r1*valc1+r2*valc2)
                ii = ii+1
            end if
!
32          continue
            call jedema()
        end do
    end if
!
! --- MENAGE
    AS_DEALLOCATE(vl=lexi)
    call jedetr('&&FOCRR2.LIR8')
    call jedetr('&&FOCRR2.INST')
!
40  continue
!
    call jedema()
end subroutine
