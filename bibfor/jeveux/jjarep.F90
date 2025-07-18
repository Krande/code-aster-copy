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
subroutine jjarep(iclas, nrmax)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=C1002
    implicit none
#include "jeveux_private.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjldyn.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jjprem.h"
#include "asterfort/jxecro.h"
#include "asterfort/jxhcod.h"
#include "asterfort/jxlibd.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: iclas, nrmax
! ----------------------------------------------------------------------
!     PERMET D'AGRANDIR UN REPERTOIRE DE NOM
!
!     IN    ICLAS  : CLASSE ASSOCIEE AU REPERTOIRE
!     IN    NRMAX  : DIMENSION DU NOUVEAU REPERTOIRE
!
! ----------------------------------------------------------------------
    integer(kind=8) :: igenr(1), itype(1), idocu(1), iorig(1), irnom(4)
    equivalence(igenr, genr), (itype, type),&
     &                 (idocu, docu), (iorig, orig), (irnom, rnom)
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic, idatin, iin, il, in, ind
    integer(kind=8) :: ipgca, iref, j, jcara, jdate, jdocu
    integer(kind=8) :: jgenr, jhcod, jiadd, jiadm, jindir, jlong
    integer(kind=8) :: jlono, jltyp, jluti, jmarq, jorig, jrnom, jtype
    integer(kind=8) :: k, khcdy, kl, kn, knody, ldynol
    integer(kind=8) :: lhcod, lloc, lon, lonoi, ltot, n, nbatt
    integer(kind=8) :: nbtot, ne
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    common/jindir/jindir(n)
!
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
! ----------------------------------------------------------------------
    character(len=32) :: clel, cle
    character(len=4) :: z
    integer(kind=8) :: jcod, kadhc, jnom, kadno, lorep, iadrs(20), kat(20)
    integer(kind=8) :: lgnom, nuti, lso(20), imq(2), iaddi(2), kdy(20)
    parameter(nbatt=12, nbtot=nbatt+3, lgnom=32)
    integer(kind=8) :: numatt(nbtot), idm(nbtot), idy(nbtot), irt
    data numatt, z/2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 16, 7, 13, 20, 17, 'INIT'/
! DEB ------------------------------------------------------------------
    ipgca = ipgc
    ipgc = -2
    irt = 0
    ic = iclas
    call jjldyn(0, -1, ltot)
!
! --- ON INTERDIT L'APPEL A JJLDYN AVEC LE PARAMETRE MODE=1 LORS DE
! --- L'ALLOCATION DYNAMIQUE  (ET LES APPELS RECURSIFS)
!
    ldynol = ldyn
    if (ldyn .eq. 1) then
        ldyn = 2
    end if
!
! --- ALLOCATION DU SEGMENT DE VALEURS POUR LE NOUVEAU REPERTOIRE
!
    lhcod = jjprem(nrmax, irt)
    if (irt .eq. 1) then
        if (ic .eq. 1) then
            call utmess('A', 'JEVEUX_64', sk=nombas(ic), si=nrmax)
        else
            call utmess('A', 'JEVEUX_65', sk=nombas(ic), si=nrmax)
        end if
    end if
    call jjalls(lhcod*lois, ic, 'V', 'I', lois, &
                'INIT', hcod, jcod, kadhc, khcdy)
    call jjecrs(kadhc, ic, 13, 0, 'E', &
                imq)
    call jjalls(nrmax*lgnom, ic, 'V', 'K', lgnom, &
                'NOINIT', irnom, jnom, kadno, knody)
    call jjecrs(kadno, ic, 7, 0, 'E', &
                imq)
    do i = 1, nrmax
        rnom(jnom-1+i) = '?'
    end do
    nuti = 0
!
! --- REMPLISSAGE DU REPERTOIRE DE NOM
!
    lorep = lhcod
    do kn = 1, nremax(ic)
        clel = rnom(jrnom(ic)+kn)
        if (clel(1:1) .eq. '?') then
            idatin = nuti+1
            goto 101
        end if
        ne = 1
        iref = jxhcod(clel, lorep)
        i = iref
5       continue
        if (hcod(jcod-1+i) .eq. 0) then
            if (nuti .ge. nrmax) then
                call utmess('F', 'JEVEUX_58')
            else
                idatin = nuti+1
                iin = i
            end if
        else
            j = hcod(jcod-1+i)
            cle = rnom(jnom-1+abs(j))
            if (cle .eq. clel) then
                call utmess('F', 'JEVEUX_59', sk=clel)
            else
                if (ne .eq. 1) in = jxhcod(clel, lorep-2)
                ne = ne+1
                i = 1+mod(i+in, lorep)
                if (ne .le. lorep) then
                    goto 5
                else
                    call utmess('F', 'JEVEUX_58')
                end if
            end if
        end if
        hcod(jcod-1+iin) = idatin
101     continue
        rnom(jnom-1+idatin) = rnom(jrnom(ic)+kn)
        nuti = nuti+1
    end do
!
! --- RECOPIE DES OBJETS SYSTEME APRES AGRANDISSEMENT
!
    nremax(ic) = nrmax
    nrhcod(ic) = lhcod
    cara(jcara(ic)) = nremax(ic)
    cara(jcara(ic)+2) = nrhcod(ic)
!
    do k = 1, nbtot
        kl = numatt(k)
        lonoi = lono(jlono(ic)+kl)*ltyp(jltyp(ic)+kl)
        if (iadd(jiadd(ic)+2*kl-1) .gt. 0) then
            call jxlibd(0, kl, ic, iadd(jiadd(ic)+2*kl-1), lonoi)
            iadd(jiadd(ic)+2*kl-1) = 0
            iadd(jiadd(ic)+2*kl) = 0
        end if
    end do
!
    kat(7) = kadno
    kdy(7) = knody
    long(jlong(ic)+7) = nrmax
    lso(7) = lhcod*lois
    kat(13) = kadhc
    kdy(13) = khcdy
    long(jlong(ic)+13) = lhcod
    lso(13) = nrmax*lgnom
!
! --- ALLOCATION DES DIFFERENTS SEGMENTS DE VALEURS :
!
    lon = 2*nremax(ic)*lois
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, imarq, iadrs(16), kat(16), kdy(16))
    call jjecrs(kat(16), ic, 16, 0, 'E', &
                imq)
    lon = nremax(ic)*len(genr(1))
    call jjalls(lon, ic, 'V', 'K', len(genr(1)), &
                z, igenr, iadrs(3), kat(3), kdy(3))
    call jjecrs(kat(3), ic, 3, 0, 'E', &
                imq)
    lon = nremax(ic)*len(type(1))
    call jjalls(lon, ic, 'V', 'K', len(type(1)), &
                z, itype, iadrs(4), kat(4), kdy(4))
    call jjecrs(kat(4), ic, 4, 0, 'E', &
                imq)
    lon = nremax(ic)*len(docu(1))
    call jjalls(lon, ic, 'V', 'K', len(docu(1)), &
                z, idocu, iadrs(5), kat(5), kdy(5))
    call jjecrs(kat(5), ic, 5, 0, 'E', &
                imq)
    lon = nremax(ic)*len(orig(1))
    call jjalls(lon, ic, 'V', 'K', len(orig(1)), &
                z, iorig, iadrs(6), kat(6), kdy(6))
    call jjecrs(kat(6), ic, 6, 0, 'E', &
                imq)
    lon = nremax(ic)*lois
    call jjalls(2*lon, ic, 'V', 'I', lois, &
                z, iadd, iadrs(2), kat(2), kdy(2))
    call jjecrs(kat(2), ic, 2, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, ltyp, iadrs(8), kat(8), kdy(8))
    call jjecrs(kat(8), ic, 8, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, long, iadrs(9), kat(9), kdy(9))
    call jjecrs(kat(9), ic, 9, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, lono, iadrs(10), kat(10), kdy(10))
    call jjecrs(kat(10), ic, 10, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, date, iadrs(11), kat(11), kdy(11))
    call jjecrs(kat(11), ic, 11, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, luti, iadrs(12), kat(12), kdy(12))
    call jjecrs(kat(12), ic, 12, 0, 'E', &
                imq)
    call jjalls(2*lon, ic, 'V', 'I', lois, &
                z, iadm, iadrs(20), kat(20), kdy(20))
    call jjecrs(kat(20), ic, 20, 0, 'E', &
                imq)
    call jjalls(lon, ic, 'V', 'I', lois, &
                z, indir, iadrs(17), kat(17), kdy(17))
    call jjecrs(kat(17), ic, 17, 0, 'E', &
                imq)
    do k = 1, nreuti(ic)
        genr(iadrs(3)-1+k) = genr(jgenr(ic)+k)
        type(iadrs(4)-1+k) = type(jtype(ic)+k)
        docu(iadrs(5)-1+k) = docu(jdocu(ic)+k)
        orig(iadrs(6)-1+k) = orig(jorig(ic)+k)
        ltyp(iadrs(8)-1+k) = ltyp(jltyp(ic)+k)
        long(iadrs(9)-1+k) = long(jlong(ic)+k)
        lono(iadrs(10)-1+k) = lono(jlono(ic)+k)
        date(iadrs(11)-1+k) = date(jdate(ic)+k)
        luti(iadrs(12)-1+k) = luti(jluti(ic)+k)
    end do
    do k = 1, 2*nreuti(ic)
        iadm(iadrs(20)-1+k) = iadm(jiadm(ic)+k)
        iadd(iadrs(2)-1+k) = iadd(jiadd(ic)+k)
        imarq(iadrs(16)-1+k) = imarq(jmarq(ic)+k)
    end do
    do i = 1, nbtot
        il = numatt(i)
        idm(i) = iadm(jiadm(ic)+2*il-1)
        idy(i) = iadm(jiadm(ic)+2*il)
    end do
    do i = 1, nbtot
        call jjlidy(idy(i), idm(i))
    end do
!
    jiadd(ic) = iadrs(2)-1
    jgenr(ic) = iadrs(3)-1
    jtype(ic) = iadrs(4)-1
    jdocu(ic) = iadrs(5)-1
    jorig(ic) = iadrs(6)-1
    jltyp(ic) = iadrs(8)-1
    jlong(ic) = iadrs(9)-1
    jlono(ic) = iadrs(10)-1
    jdate(ic) = iadrs(11)-1
    jluti(ic) = iadrs(12)-1
    jmarq(ic) = iadrs(16)-1
    kmarq(ic) = kat(16)
    jiadm(ic) = iadrs(20)-1
    kiadm(ic) = kat(20)
    jhcod(ic) = jcod-1
    jrnom(ic) = jnom-1
    jindir(ic) = iadrs(17)-1
!
    do i = 1, nbtot
        il = numatt(i)
        if (il .eq. 13) then
            long(jlong(ic)+il) = lhcod
            lono(jlono(ic)+il) = lhcod
        else if (il .eq. 16 .or. il .eq. 2 .or. il .eq. 20) then
            long(jlong(ic)+il) = 2*nrmax
            lono(jlono(ic)+il) = 2*nrmax
        else
            long(jlong(ic)+il) = nrmax
            lono(jlono(ic)+il) = nrmax
        end if
        lloc = lono(jlono(ic)+il)*ltyp(jltyp(ic)+il)
        if (mod(lloc, lois) .ne. 0) then
            lono(jlono(ic)+il) = ((1+lloc/lois)*lois)/ltyp(jltyp(ic)+il)
        end if
        lso(il) = lono(jlono(ic)+il)*ltyp(jltyp(ic)+il)
        iadm(jiadm(ic)+2*il-1) = kat(il)
        iadm(jiadm(ic)+2*il) = kdy(il)
    end do
!
    do i = 1, nbtot
        il = numatt(i)
        iaddi(1) = 0
        iaddi(2) = 0
        call jxecro(ic, kat(il), iaddi, lso(il), 0, &
                    il)
        iadd(jiadd(ic)+2*il-1) = iaddi(1)
        iadd(jiadd(ic)+2*il) = iaddi(2)
    end do
    cara(jcara(ic)+6) = iadd(jiadd(ic)+2*2-1)
    cara(jcara(ic)+7) = iadd(jiadd(ic)+2*2)
    do ind = 1, nremax(ic)
        indir(jindir(ic)+ind) = ind
    end do
    ldyn = ldynol
    ipgc = ipgca
! FIN ------------------------------------------------------------------
end subroutine
