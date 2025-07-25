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
subroutine jeinif(sti, sto, nomf, clas, nrep, &
                  nbloc, lbloc)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=C1002
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/gtopti.h"
#include "asterc/gtoptk.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjprem.h"
#include "asterfort/jxecro.h"
#include "asterfort/jxlibd.h"
#include "asterfort/jxlir1.h"
#include "asterfort/jxliro.h"
#include "asterfort/jxouvr.h"
#include "asterfort/lxmins.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nrep, nbloc, lbloc
    character(len=*) :: sti, sto, nomf, clas
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR D'OUVERTURE D'UNE CLASSE
!
! IN  STI    : STATUT EN DEBUT DE TRAVAIL ('DEBUT','POURSUIT')
! IN  STO    : STATUT EN FIN DE TRAVAIL ('SAUVE','DETRUIT')
! IN  NOMF   : NOM LOCALE DE LA BASE
! IN  CLAS   : NOM LOCALE DE LA CLASSE
! IN  NREP   : LONGUEUR DU REPERTOIRE
! IN  NBLOC  : NOMBRE D'ENREGISTREMMENTS DU FICHIER D'ACCES DIRECT
!              SI NBLOC = 0,  ON LE DETERMINE A PARTIR DE MFIC
! IN  LBLOC  : LONGUEUR DES ENREGISTREMMENTS DU FICHIER D'ACCES DIRECT
!              CETTE LONGUEUR EST DONNEE EN KILO (1024) MOT (ENTIER)
!
! ----------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadrs, ic, icre, ipgca
    integer(kind=8) :: iret, jcara, jdate, jdocu, jgenr, jhcod
    integer(kind=8) :: jiacce, jiadd, jiadm, jindir, jlong, jlono, jltyp
    integer(kind=8) :: jluti, jmarq, jorig, jrnom, jtype, jusadi, k
    integer(kind=8) :: l, lcarao, ldynol, lloc, lmarq, lon, lon1
    integer(kind=8) :: lon2, n, nbacce, nbext, nbgros, nblim, nblma1
    integer(kind=8) :: nblma2, nbloco, nbpeti, nrepo, mfic_prev
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    common/jiacce/jiacce(n), nbacce(2*n)
    common/jusadi/jusadi(n)
    common/jindir/jindir(n)
    common/inbdet/nblim(n), nbgros(n), nbpeti(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
!
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    integer(kind=8) :: idn, iext, nbenrg
    common/iextje/idn(n), iext(n), nbenrg(n)
    integer(kind=8) :: nbcla
    common/nficje/nbcla
! ----------------------------------------------------------------------
    integer(kind=8) :: igenr(1), itype(1), idocu(1), iorig(1), irnom(4)
    equivalence(igenr, genr), (itype, type),&
     &                 (idocu, docu), (iorig, orig), (irnom, rnom)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: lfic, mfic
    common/fenvje/lfic(n), mfic
!
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
! ----------------------------------------------------------------------
    character(len=1) :: kclas
    character(len=4) :: z
    parameter(z='INIT')
    character(len=8) :: knom, knomf, kstin, kstou, cversb, cversu
    character(len=24) :: valk(3)
    integer(kind=8) :: ncar, itlec(1), itecr(1), iadadd(2), lgbl
    integer(kind=8) :: vali(7), irt, ind, iesup
    parameter(ncar=13)
! ----------------------------------------------------------------------
    aster_logical :: lenrg
    integer(kind=8) :: lidbas, lideff
    parameter(lidbas=20, lideff=15)
    character(len=8) :: cidbas(lidbas)
    integer(kind=8) :: kat(lidbas), lso(lidbas), kdy(lidbas)
    data cidbas/'$$CARA  ', '$$IADD  ', '$$GENR  ', '$$TYPE  ',&
     &               '$$DOCU  ', '$$ORIG  ', '$$RNOM  ', '$$LTYP  ',&
     &               '$$LONG  ', '$$LONO  ', '$$DATE  ', '$$LUTI  ',&
     &               '$$HCOD  ', '$$USADI ', '$$ACCE  ', '$$MARQ  ',&
     &               '$$INDI  ', '$$TLEC  ', '$$TECR  ', '$$IADM  '/
! DEB ------------------------------------------------------------------
    ipgca = ipgc
    ipgc = -2
    irt = 0
!
    kclas = clas
    kstin = sti
    kstou = sto
    knom = nomf
    call lxmins(nomf)
    knomf = nomf
!
    ASSERT(knomf .ne. '        ' .and. len(nomf) .le. 8)
    ASSERT(kclas .ne. ' ')
    ASSERT(index(classe, kclas) .eq. 0)
!
    ASSERT(kstin .eq. 'DEBUT   ' .or. kstin .eq. 'POURSUIT')
    ASSERT(kstou .eq. 'SAUVE   ' .or. kstou .eq. 'DETRUIT ' .or. kstou .eq. 'LIBERE ')
    ASSERT(nrep .gt. 0)
    ASSERT(lbloc .gt. 0)
!
    ic = index(classe, ' ')
    ASSERT(ic .gt. 0)
    nomfic(ic) = knomf
    nombas(ic) = knom
    kstini(ic) = kstin
    kstout(ic) = kstou
    classe(ic:ic) = kclas
    nbcla = index(classe, '$')-1
    if (nbcla .eq. -1) nbcla = n
!
    iclas = ic
    nbgros(ic) = 0
    nbpeti(ic) = 0
    nomuti = ' '
    nomos = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
    nomco = '$$$$$$$$$$$$$$$$$$$$$$$$'
    nomoc = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!
! --- ON INTERDIT L'APPEL A JJLDYN LORS DE L'ALLOCATION
! --- DYNAMIQUE  (ET LES APPELS RECURSIFS)
!
    ldynol = ldyn
    if (ldyn .eq. 1) then
        ldyn = 2
    end if
!
    if (kstin .eq. 'DEBUT   ') then
        nremax(ic) = nrep
        nreuti(ic) = 0
        nrhcod(ic) = jjprem(nremax(ic), irt)
        if (irt .eq. 1) then
            if (ic .eq. 1) then
                call utmess('A', 'JEVEUX_64', sk=nombas(ic), si=nremax(ic))
            else
                call utmess('A', 'JEVEUX_65', sk=nombas(ic), si=nremax(ic))
            end if
        end if
        nbluti(ic) = 0
        longbl(ic) = lbloc
        if (nbloc .eq. 0) then
            nblmax(ic) = mfic/(longbl(ic)*lois)
        else
            nblmax(ic) = min(nbloc, mfic/(longbl(ic)*lois))
        end if
        nblim(ic) = 500
!
        lmarq = 2*nrep*lois
        call jjalls(lmarq, ic, 'V', 'I', lois, &
                    z, imarq, iadrs, kmarq(ic), kdy(16))
        kat(16) = kmarq(ic)
        jmarq(ic) = iadrs-1
        call jjecrs(kat(16), ic, 16, 0, 'E', &
                    imarq(jmarq(ic)+2*16-1))
!
        lcarao = ncar*lois
        call jjalls(lcarao, ic, 'V', 'I', lois, &
                    z, cara, iadrs, kat(1), kdy(1))
        jcara(ic) = iadrs
        call jjecrs(kat(1), ic, 1, 0, 'E', &
                    imarq(jmarq(ic)+1))
!
        nbenrg(ic) = min(lfic(ic)/(longbl(ic)*lois), nblmax(ic))
!
        nbloco = nblmax(ic)*lois
        call jjalls(nbloco, ic, 'V', 'I', lois, &
                    z, iacce, iadrs, kat(15), kdy(15))
        jiacce(ic) = iadrs-1
        call jjecrs(kat(15), ic, 15, 0, 'E', &
                    imarq(jmarq(ic)+2*15-1))
!
        lgbl = nremax(ic)*lois
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, indir, iadrs, kat(17), kdy(17))
        jindir(ic) = iadrs-1
        call jjecrs(kat(17), ic, 17, 0, 'E', &
                    imarq(jmarq(ic)+2*17-1))
        do ind = 1, nremax(ic)
            indir(jindir(ic)+ind) = ind
        end do
!
        lgbl = 1024*longbl(ic)*lois
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, itlec, iadrs, kitlec(ic), kdy(18))
        kat(18) = kitlec(ic)
        kitlec(ic) = (kitlec(ic)-1)*lois
        call jjecrs(kat(18), ic, 18, 0, 'E', &
                    imarq(jmarq(ic)+2*18-1))
!
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, itecr, iadrs, kitecr(ic), kdy(19))
        kat(19) = kitecr(ic)
        kitecr(ic) = (kitecr(ic)-1)*lois
        call jjecrs(kat(19), ic, 19, 0, 'E', &
                    imarq(jmarq(ic)+2*19-1))
!
        nrepo = nrep*lois
        call jjalls(2*nrepo, ic, 'V', 'I', lois, &
                    z, iadm, iadrs, kiadm(ic), kdy(20))
        kat(20) = kiadm(ic)
        jiadm(ic) = iadrs-1
        call jjecrs(kat(20), ic, 20, 0, 'E', &
                    imarq(jmarq(ic)+2*20-1))
!
! ----- OPEN DU FICHIER
!
        call jxouvr(ic, 1)
        iext(ic) = 1
!
! ----- ECRITURE DANS L'OBJET CARA
! ----- NREMAX NREUTI NRHCOD NBLMAX NBLUTI
!
        cara(jcara(ic)) = nremax(ic)
        cara(jcara(ic)+1) = nreuti(ic)
        cara(jcara(ic)+2) = nrhcod(ic)
        cara(jcara(ic)+3) = nblmax(ic)
        cara(jcara(ic)+4) = nbluti(ic)
        cara(jcara(ic)+5) = longbl(ic)
        call gtopti('versMAJ', cara(jcara(ic)+8), iret)
        call gtopti('versMIN', cara(jcara(ic)+9), iret)
        call gtopti('versSUB', cara(jcara(ic)+10), iret)
        cara(jcara(ic)+11) = lfic(ic)
        lon = 2*nremax(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, iadd, iadrs, kat(2), kdy(2))
        jiadd(ic) = iadrs-1
        call jjecrs(kat(2), ic, 2, 0, 'E', &
                    imarq(jmarq(ic)+2*2-1))
        lon = nremax(ic)*len(genr(1))
        call jjalls(lon, ic, 'V', 'K', len(genr(1)), &
                    z, igenr, iadrs, kat(3), kdy(3))
        jgenr(ic) = iadrs-1
        call jjecrs(kat(3), ic, 3, 0, 'E', &
                    imarq(jmarq(ic)+2*3-1))
        lon = nremax(ic)*len(type(1))
        call jjalls(lon, ic, 'V', 'K', len(type(1)), &
                    z, itype, iadrs, kat(4), kdy(4))
        jtype(ic) = iadrs-1
        call jjecrs(kat(4), ic, 4, 0, 'E', &
                    imarq(jmarq(ic)+2*4-1))
        lon = nremax(ic)*len(docu(1))
        call jjalls(lon, ic, 'V', 'K', len(docu(1)), &
                    z, idocu, iadrs, kat(5), kdy(5))
        jdocu(ic) = iadrs-1
        call jjecrs(kat(5), ic, 5, 0, 'E', &
                    imarq(jmarq(ic)+2*5-1))
        lon = nremax(ic)*len(orig(1))
        call jjalls(lon, ic, 'V', 'K', len(orig(1)), &
                    z, iorig, iadrs, kat(6), kdy(6))
        jorig(ic) = iadrs-1
        call jjecrs(kat(6), ic, 6, 0, 'E', &
                    imarq(jmarq(ic)+2*6-1))
        lon = nremax(ic)*len(rnom(1))
        call jjalls(lon, ic, 'V', 'K', len(rnom(1)), &
                    z, irnom, iadrs, kat(7), kdy(7))
        jrnom(ic) = iadrs-1
        call jjecrs(kat(7), ic, 7, 0, 'E', &
                    imarq(jmarq(ic)+2*7-1))
        do ind = 1, nremax(ic)
            rnom(jrnom(ic)+ind) = '?'
        end do
        lon = nremax(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, ltyp, iadrs, kat(8), kdy(8))
        jltyp(ic) = iadrs-1
        call jjecrs(kat(8), ic, 8, 0, 'E', &
                    imarq(jmarq(ic)+2*8-1))
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, long, iadrs, kat(9), kdy(9))
        jlong(ic) = iadrs-1
        call jjecrs(kat(9), ic, 9, 0, 'E', &
                    imarq(jmarq(ic)+2*9-1))
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, lono, iadrs, kat(10), kdy(10))
        jlono(ic) = iadrs-1
        call jjecrs(kat(10), ic, 10, 0, 'E', &
                    imarq(jmarq(ic)+2*10-1))
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, date, iadrs, kat(11), kdy(11))
        jdate(ic) = iadrs-1
        call jjecrs(kat(11), ic, 11, 0, 'E', &
                    imarq(jmarq(ic)+2*11-1))
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, luti, iadrs, kat(12), kdy(12))
        jluti(ic) = iadrs-1
        call jjecrs(kat(12), ic, 12, 0, 'E', &
                    imarq(jmarq(ic)+2*12-1))
        lon = nrhcod(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, hcod, iadrs, kat(13), kdy(13))
        jhcod(ic) = iadrs-1
        call jjecrs(kat(13), ic, 13, 0, 'E', &
                    imarq(jmarq(ic)+2*13-1))
        lon = 3*nblmax(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, iusadi, iadrs, kat(14), kdy(14))
        do l = 1, nblmax(ic)
            iusadi(iadrs+(3*l-2)-1) = -1
            iusadi(iadrs+(3*l-1)-1) = -1
            iusadi(iadrs+(3*l)-1) = 0
        end do
        jusadi(ic) = iadrs-1
        call jjecrs(kat(14), ic, 14, 0, 'E', &
                    imarq(jmarq(ic)+2*14-1))
        icre = 1
        ltyp(jltyp(ic)+3) = len(genr(1))
        ltyp(jltyp(ic)+4) = len(type(1))
        ltyp(jltyp(ic)+5) = len(docu(1))
        ltyp(jltyp(ic)+6) = len(orig(1))
        ltyp(jltyp(ic)+7) = len(rnom(1))
        do i = 1, lidbas
            nomuti = '________'//nombas(ic)//'________'//cidbas(i)
            call jjcren(nomuti, icre, iret)
            genr(jgenr(ic)+i) = 'V'
            if ((i .ge. 3 .and. i .le. 7)) then
                type(jtype(ic)+i) = 'K'
            else
                type(jtype(ic)+i) = 'I'
                ltyp(jltyp(ic)+i) = lois
            end if
            if (i .eq. 1) then
                long(jlong(ic)+i) = ncar
                lono(jlono(ic)+i) = ncar
                iadd(jiadd(ic)+2*i-1) = 0
                iadd(jiadd(ic)+2*i) = 0
                call jxecro(ic, kat(1), iadd(jiadd(ic)+2*i-1), lono(jlono(ic)+i)*lois, 0, &
                            1)
            else if (i .eq. 2 .or. i .eq. 16 .or. i .eq. 20) then
                long(jlong(ic)+i) = 2*nremax(ic)
                lono(jlono(ic)+i) = 2*nremax(ic)
                lso(i) = 2*nremax(ic)*lois
            else if (i .eq. 13) then
                long(jlong(ic)+i) = nrhcod(ic)
                lono(jlono(ic)+i) = nrhcod(ic)
                lso(i) = nrhcod(ic)*ltyp(jltyp(ic)+i)
            else if (i .eq. 14) then
                long(jlong(ic)+i) = 3*nblmax(ic)
                lono(jlono(ic)+i) = 3*nblmax(ic)
                lso(i) = 3*nblmax(ic)*ltyp(jltyp(ic)+i)
            else if (i .eq. 15) then
                long(jlong(ic)+i) = nblmax(ic)
                lono(jlono(ic)+i) = nblmax(ic)
                lso(i) = nblmax(ic)*ltyp(jltyp(ic)+i)
            else if (i .eq. 18 .or. i .eq. 19) then
                long(jlong(ic)+i) = 1024*longbl(ic)
                lono(jlono(ic)+i) = 1024*longbl(ic)
                lso(i) = 1024*longbl(ic)*ltyp(jltyp(ic)+i)
            else
                long(jlong(ic)+i) = nremax(ic)
                lono(jlono(ic)+i) = nremax(ic)
                lloc = lono(jlono(ic)+i)*ltyp(jltyp(ic)+i)
                if (mod(lloc, lois) .ne. 0) then
                    lono(jlono(ic)+i) = ((1+lloc/lois)*lois)/ltyp(jltyp(ic)+i)
                end if
                lso(i) = lono(jlono(ic)+i)*ltyp(jltyp(ic)+i)
            end if
            iadm(jiadm(ic)+2*i-1) = kat(i)
            iadm(jiadm(ic)+2*i) = kdy(i)
        end do
!
        do i = 2, lideff
            iadd(jiadd(ic)+2*i-1) = 0
            iadd(jiadd(ic)+2*i) = 0
            call jxecro(ic, kat(i), iadd(jiadd(ic)+2*i-1), lso(i), 0, &
                        i)
        end do
        cara(jcara(ic)+6) = iadd(jiadd(ic)+2*2-1)
        cara(jcara(ic)+7) = iadd(jiadd(ic)+2*2)
    else
!
! ----- OPEN FICHIER
! ----- LECTURE DANS LE PREMIER BLOC DU FICHIER ET FERMETURE
!
        lcarao = ncar*lois
        call jjalls(lcarao, ic, 'V', 'I', lois, &
                    z, cara, iadrs, kat(1), kdy(1))
        jcara(ic) = iadrs
!
!  ---- L'ECRITURE DU STATUT ET DE L'ETAT SERA DE NOUVEAU EFFECTUEE
! ----- LORSQUE LES DIMENSIONS AURONT ETE RELUES
!
        call jjecrs(kat(1), ic, 1, 0, 'E', &
                    vali)
!
        call jxlir1(ic, cara(jcara(ic)))
        cversb = '  .  .  '
        call codent(cara(jcara(ic)+8), 'D ', cversb(1:2))
        call codent(cara(jcara(ic)+9), 'D0', cversb(4:5))
        call codent(cara(jcara(ic)+10), 'D0', cversb(7:8))
        call gtoptk('versionD0', cversu, iret)
        nremax(ic) = cara(jcara(ic))
        nreuti(ic) = cara(jcara(ic)+1)
        nrhcod(ic) = cara(jcara(ic)+2)
        nblmax(ic) = cara(jcara(ic)+3)
        nbluti(ic) = cara(jcara(ic)+4)
        longbl(ic) = cara(jcara(ic)+5)
        iadadd(1) = cara(jcara(ic)+6)
        iadadd(2) = cara(jcara(ic)+7)
        lfic(ic) = cara(jcara(ic)+11)
        nbenrg(ic) = cara(jcara(ic)+12)
        if (cversu .ne. cversb) then
            valk(1) = nombas(ic)
            valk(2) = cversb
            valk(3) = cversu
            call utmess('A', 'JEVEUX_08', nk=3, valk=valk)
        end if
!
!       calcul de la taille maximale précédemment utilisée
        mfic_prev = longbl(ic)*nblmax(ic)*lois
        if (mfic_prev .ne. mfic .and. iext(ic) .gt. 0) then
            mfic = mfic_prev
            call utmess("I", "JEVEUX1_79", sr=1.0d0*mfic/1024/1024)
        end if
!
        if (nbloc .eq. 0) then
            nblma2 = mfic/(longbl(ic)*lois)
        else
            nblma2 = min(nbloc, mfic/(longbl(ic)*lois))
        end if
!
! ---- LORSQUE LE NOMBRE D'ENREGISTREMENTS MAXIMUM EST MODIFIE
!
        nblma1 = nblmax(ic)
        if (nblmax(ic) .ge. nblma2) then
            lenrg = .false.
            nblma2 = nblmax(ic)
        else
            vali(1) = nblmax(ic)
            vali(2) = nblma2
            valk(1) = nombas(ic)
            call utmess('I', 'JEVEUX_36', sk=valk(1), ni=2, vali=vali)
            lenrg = .true.
        end if
!
        valk(1) = nombas(ic)
        valk(2) = cversb
        vali(1) = nbluti(ic)
        vali(2) = nblmax(ic)
        vali(3) = 1024*longbl(ic)*lois
        vali(4) = nreuti(ic)
        vali(5) = nremax(ic)
        vali(6) = (nreuti(ic)*100)/nremax(ic)
        vali(7) = nbenrg(ic)
!
        call utmess('I', 'JEVEUX_21', nk=2, valk=valk, ni=7, &
                    vali=vali)
!
        nblmax(ic) = nblma2
!
        lmarq = 2*nremax(ic)*lois
        call jjalls(lmarq, ic, 'V', 'I', lois, &
                    z, imarq, iadrs, kmarq(ic), kdy(16))
        kat(16) = kmarq(ic)
        jmarq(ic) = iadrs-1
        call jjecrs(kat(16), ic, 16, 0, 'E', &
                    imarq(jmarq(ic)+2*16-1))
!
        call jjecrs(kat(1), ic, 1, 0, 'E', &
                    imarq(jmarq(ic)+2*1-1))
!
        if (nbenrg(ic) .eq. 0) then
            nbenrg(ic) = min(lfic(ic)/(longbl(ic)*lois), nblma2)
        end if
!
! ----- NOUVEL OPEN DE LA BASE
        iesup = 1
        if (mod(nbluti(ic), nbenrg(ic)) .eq. 0) iesup = 0
        nbext = (nbluti(ic)/nbenrg(ic))+iesup
        do k = 0, nbext-1
            call jxouvr(ic, k+1)
        end do
        iext(ic) = nbext
!
        lgbl = nremax(ic)*lois
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, indir, iadrs, kat(17), kdy(17))
        jindir(ic) = iadrs-1
        call jjecrs(kat(17), ic, 17, 0, 'E', &
                    imarq(jmarq(ic)+2*17-1))
        do ind = 1, nremax(ic)
            indir(jindir(ic)+ind) = ind
        end do
!
        lgbl = 1024*longbl(ic)*lois
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, itlec, iadrs, kitlec(ic), kdy(18))
        kat(18) = kitlec(ic)
        kitlec(ic) = (kitlec(ic)-1)*lois
        call jjecrs(kat(18), ic, 18, 0, 'E', &
                    imarq(jmarq(ic)+2*18-1))
        call jjalls(lgbl, ic, 'V', 'I', lois, &
                    z, itecr, iadrs, kitecr(ic), kdy(19))
        kat(19) = kitecr(ic)
        kitecr(ic) = (kitecr(ic)-1)*lois
        call jjecrs(kat(19), ic, 19, 0, 'E', &
                    imarq(jmarq(ic)+2*19-1))
        lon = nremax(ic)*lois
        call jjalls(2*lon, ic, 'V', 'I', lois, &
                    z, iadm, iadrs, kiadm(ic), kdy(20))
        kat(20) = kiadm(ic)
        jiadm(ic) = iadrs-1
        call jjecrs(kat(20), ic, 20, 0, 'E', &
                    imarq(jmarq(ic)+2*20-1))
!
        lon2 = nblma2*lois
        call jjalls(lon2, ic, 'V', 'I', lois, &
                    z, iacce, iadrs, kat(15), kdy(15))
        jiacce(ic) = iadrs-1
        call jjecrs(kat(15), ic, 15, 0, 'E', &
                    imarq(jmarq(ic)+2*15-1))
!
        call jjalls(2*lon, ic, 'V', 'I', lois, &
                    z, iadd, iadrs, kat(2), kdy(2))
        jiadd(ic) = iadrs-1
        call jjecrs(kat(2), ic, 2, 0, 'E', &
                    imarq(jmarq(ic)+2*2-1))
!
        call jxliro(ic, kat(2), iadadd, 2*lon)
!
        lon2 = 3*nblma2*lois
        call jjalls(lon2, ic, 'V', 'I', lois, &
                    z, iusadi, iadrs, kat(14), kdy(14))
        do l = 1, nblma2
            iusadi(iadrs+(3*l-2)-1) = -1
            iusadi(iadrs+(3*l-1)-1) = -1
            iusadi(iadrs+(3*l)-1) = 0
        end do
        jusadi(ic) = iadrs-1
        call jjecrs(kat(14), ic, 14, 0, 'E', &
                    imarq(jmarq(ic)+2*14-1))
        lon1 = 3*nblma1*lois
        call jxliro(ic, kat(14), iadd(jiadd(ic)+2*14-1), lon1)
        if (lenrg) then
            call jxlibd(0, 14, ic, iadd(jiadd(ic)+2*14-1), lon1)
            iadd(jiadd(ic)+2*14-1) = 0
            iadd(jiadd(ic)+2*14) = 0
            call jxecro(ic, kat(14), iadd(jiadd(ic)+2*14-1), lon2, 0, &
                        14)
        end if
!
        lon1 = nblma1*lois
        call jxliro(ic, kat(15), iadd(jiadd(ic)+2*15-1), lon1)
        if (lenrg) then
            call jxlibd(0, 15, ic, iadd(jiadd(ic)+2*15-1), lon1)
            iadd(jiadd(ic)+2*15-1) = 0
            iadd(jiadd(ic)+2*15) = 0
            call jxecro(ic, kat(15), iadd(jiadd(ic)+2*15-1), lon2, 0, &
                        15)
        end if
!
        lon = nremax(ic)*len(genr(1))
        call jjalls(lon, ic, 'V', 'K', len(genr(1)), &
                    z, igenr, iadrs, kat(3), kdy(3))
        jgenr(ic) = iadrs-1
        call jjecrs(kat(3), ic, 3, 0, 'E', &
                    imarq(jmarq(ic)+2*3-1))
        call jxliro(ic, kat(3), iadd(jiadd(ic)+2*3-1), lon)
!
        lon = nremax(ic)*len(type(1))
        call jjalls(lon, ic, 'V', 'K', len(type(1)), &
                    z, itype, iadrs, kat(4), kdy(4))
        jtype(ic) = iadrs-1
        call jjecrs(kat(4), ic, 4, 0, 'E', &
                    imarq(jmarq(ic)+2*4-1))
        call jxliro(ic, kat(4), iadd(jiadd(ic)+2*4-1), lon)
        lon = nremax(ic)*len(docu(1))
        call jjalls(lon, ic, 'V', 'K', len(docu(1)), &
                    z, idocu, iadrs, kat(5), kdy(5))
        jdocu(ic) = iadrs-1
        call jjecrs(kat(5), ic, 5, 0, 'E', &
                    imarq(jmarq(ic)+2*5-1))
        call jxliro(ic, kat(5), iadd(jiadd(ic)+2*5-1), lon)
        lon = nremax(ic)*len(orig(1))
        call jjalls(lon, ic, 'V', 'K', len(orig(1)), &
                    z, iorig, iadrs, kat(6), kdy(6))
        jorig(ic) = iadrs-1
        call jjecrs(kat(6), ic, 6, 0, 'E', &
                    imarq(jmarq(ic)+2*6-1))
        call jxliro(ic, kat(6), iadd(jiadd(ic)+2*6-1), lon)
        lon = nremax(ic)*len(rnom(1))
        call jjalls(lon, ic, 'V', 'K', len(rnom(1)), &
                    z, irnom, iadrs, kat(7), kdy(7))
        jrnom(ic) = iadrs-1
        call jjecrs(kat(7), ic, 7, 0, 'E', &
                    imarq(jmarq(ic)+2*7-1))
        call jxliro(ic, kat(7), iadd(jiadd(ic)+2*7-1), lon)
        lon = nremax(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, ltyp, iadrs, kat(8), kdy(8))
        jltyp(ic) = iadrs-1
        call jjecrs(kat(8), ic, 8, 0, 'E', &
                    imarq(jmarq(ic)+2*8-1))
        call jxliro(ic, kat(8), iadd(jiadd(ic)+2*8-1), lon)
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, long, iadrs, kat(9), kdy(9))
        jlong(ic) = iadrs-1
        call jjecrs(kat(9), ic, 9, 0, 'E', &
                    imarq(jmarq(ic)+2*9-1))
        call jxliro(ic, kat(9), iadd(jiadd(ic)+2*9-1), lon)
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, lono, iadrs, kat(10), kdy(10))
        jlono(ic) = iadrs-1
        call jjecrs(kat(10), ic, 10, 0, 'E', &
                    imarq(jmarq(ic)+2*10-1))
        call jxliro(ic, kat(10), iadd(jiadd(ic)+2*10-1), lon)
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, date, iadrs, kat(11), kdy(11))
        jdate(ic) = iadrs-1
        call jjecrs(kat(11), ic, 11, 0, 'E', &
                    imarq(jmarq(ic)+2*11-1))
        call jxliro(ic, kat(11), iadd(jiadd(ic)+2*11-1), lon)
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, luti, iadrs, kat(12), kdy(12))
        jluti(ic) = iadrs-1
        call jjecrs(kat(12), ic, 12, 0, 'E', &
                    imarq(jmarq(ic)+2*12-1))
        call jxliro(ic, kat(12), iadd(jiadd(ic)+2*12-1), lon)
        lon = nrhcod(ic)*lois
        call jjalls(lon, ic, 'V', 'I', lois, &
                    z, hcod, iadrs, kat(13), kdy(13))
        jhcod(ic) = iadrs-1
        call jjecrs(kat(13), ic, 13, 0, 'E', &
                    imarq(jmarq(ic)+2*13-1))
        call jxliro(ic, kat(13), iadd(jiadd(ic)+2*13-1), lon)
        do i = 1, lidbas
            iadm(jiadm(ic)+2*i-1) = kat(i)
            iadm(jiadm(ic)+2*i) = kdy(i)
        end do
        if (lenrg) then
            long(jlong(ic)+15) = nblma2
            lono(jlono(ic)+15) = nblma2
            long(jlong(ic)+14) = 3*nblma2
            lono(jlono(ic)+14) = 3*nblma2
            lon2 = lono(jlono(ic)+9)*ltyp(jltyp(ic)+9)
            call jxecro(ic, kat(9), iadd(jiadd(ic)+2*9-1), lon2, 0, &
                        9)
            lon2 = lono(jlono(ic)+10)*ltyp(jltyp(ic)+10)
            call jxecro(ic, kat(10), iadd(jiadd(ic)+2*10-1), lon2, 0, &
                        10)
            call jxecro(ic, kat(1), iadd(jiadd(ic)+2*1-1), lono(jlono(ic)+1)*lois, 0, &
                        1)
        end if
    end if
!
    ldyn = ldynol
    ipgc = ipgca
! FIN ------------------------------------------------------------------
end subroutine
