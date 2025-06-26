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
subroutine jelibf(cond, clas, info)
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR PERMETTANT DE LIBERER TOUS LES OBJETS D'UNE BASE
!         ET DE FERMER LES FICHIERS ASSOCIES
! IN : CLAS NOM DE LA CLASSE ASSOCIEE ('G','V', ...)
! IN : COND TYPE DE FERMETURE
!      = 'SAUVE'
!      = 'ERREUR'
!      = 'DETRUIT'
!      = 'LIBERE'
! IN : INFO = 0 MODE SILENCIEUX
!      INFO = 1 MODE BAVARD
!
! ----------------------------------------------------------------------
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/rmfile.h"
#include "asterfort/assert.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jxecrb.h"
#include "asterfort/jxecro.h"
#include "asterfort/jxferm.h"
#include "asterfort/lxmins.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: cond
    character(len=*), intent(in) :: clas
    integer(kind=8), intent(in) :: info
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad2, iadacc, iadacy, iadadi, iadady, infrm
    integer(kind=8) :: iadmi, iadyn, ibacol, ic, idb, iret
    integer(kind=8) :: jcara, jdate, jdocu, jgenr, jhcod, jiacce, jiadd
    integer(kind=8) :: jiadm, jindir, jlong, jlono, jltyp, jluti, jmarq
    integer(kind=8) :: jorig, jrnom, jtype, jusadi, k16, k17, k17i
    integer(kind=8) :: k18, k18i, k19, k19i, k2, k20, lacc
    integer(kind=8) :: ladi, ladi2, n, nbacce, iadad2
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    aster_logical :: litlec
    common/lficje/litlec(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n), dn2(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    integer(kind=8) :: idn, iext, nbenrg
    common/iextje/idn(n), iext(n), nbenrg(n)
!
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
!
    common/jiacce/jiacce(n), nbacce(2*n)
    common/jusadi/jusadi(n)
    common/jindir/jindir(n)
    integer(kind=8) :: lundef, idebug
    common/undfje/lundef, idebug
! ----------------------------------------------------------------------
    integer(kind=8) :: lidbas, lideff
    parameter(lidbas=20, lideff=15)
    character(len=1) :: kclas
    character(len=8) :: kcond, valk(1), nom
    character(len=32) :: nomcar
    integer(kind=8) :: iadcar, iaddac(2), lgbl, vali(9), iadcdy, iaddad(2)
    real(kind=8) :: valr(2)
    aster_logical :: bool
! DEB ------------------------------------------------------------------
    kcond = cond
    kclas = clas
    ic = index(classe, kclas)
    iclaos = ic
    ASSERT(ic .ne. 0)
    bool = kcond .eq. '        ' .or. kcond .eq. 'SAUVE   ' .or. kcond .eq. 'ERREUR  ' .or. kcond &
           .eq. 'DETRUIT ' .or. kcond .eq. 'LIBERE  '
    ASSERT(bool)
    if (kcond .eq. '        ') then
        kcond = kstout(ic)
    else if (kcond .eq. 'ERREUR  ') then
        kcond = 'SAUVE'
    end if
!
!     ----------- DECHARGER TOUTES LES COLLECTIONS -----------------
    do i = lidbas+1, nreuti(ic)
        if (genr(jgenr(ic)+i) .eq. 'X') then
            iclaco = ic
            ibacol = iadm(jiadm(ic)+2*i-1)
            if (ibacol .ne. 0) then
                idatco = i
                nomco = rnom(jrnom(ic)+i) (1:24)
                call jjlide('JELIBF', nomco, 2)
            end if
        end if
    end do
!     -- DECHARGER TOUS LES OBJETS SIMPLES Y COMPRIS AVEC $$ -------
    do i = lidbas+1, nreuti(ic)
        iadmi = iadm(jiadm(ic)+2*i-1)
        if (iadmi .ne. 0) then
            idatos = i
            nomos = rnom(jrnom(ic)+i)
            call jjlide('JELIBF', nomos, 1)
        end if
    end do
!     ----------- DECHARGER TOUS LES OBJETS SYSTEME ----------------
    iad2 = iadm(jiadm(ic)+2*2-1)
    iadcar = iadm(jiadm(ic)+2*1-1)
    iadcdy = iadm(jiadm(ic)+2*1)
    iadacc = iadm(jiadm(ic)+2*lideff-1)
    iadacy = iadm(jiadm(ic)+2*lideff)
    k2 = iadm(jiadm(ic)+2*2)
    k16 = iadm(jiadm(ic)+2*16)
    k17 = iadm(jiadm(ic)+2*17)
    k17i = iadm(jiadm(ic)+2*17-1)
    k18 = iadm(jiadm(ic)+2*18)
    k18i = iadm(jiadm(ic)+2*18-1)
    k19 = iadm(jiadm(ic)+2*19)
    k19i = iadm(jiadm(ic)+2*19-1)
    k20 = iadm(jiadm(ic)+2*20)
!
    if (kcond .ne. 'LIBERE  ') then
!
        nomcar = rnom(jrnom(ic)+1)
        iadcar = iadm(jiadm(ic)+2*1-1)
        lacc = lono(jlono(ic)+lideff)*ltyp(jltyp(ic)+lideff)
        iaddac(1) = iadd(jiadd(ic)+2*lideff-1)
        iaddac(2) = iadd(jiadd(ic)+2*lideff)
!
        iadadi = iadm(jiadm(ic)+2*(lideff-1)-1)
        iadady = iadm(jiadm(ic)+2*(lideff-1))
        ladi = lono(jlono(ic)+lideff-1)*ltyp(jltyp(ic)+lideff-1)
        iaddad(1) = iadd(jiadd(ic)+2*(lideff-1)-1)
        iaddad(2) = iadd(jiadd(ic)+2*(lideff-1))
!
! --- ON ECRIT LES OBJETS SYSTEME ($$USADI idos=14 ET $$ACCE idos=15) DONT LES LONGUEURS ONT PU ETRE
! --- MODIFIEES LORS DE L'APPEL A jjagod (REDIMENSIONNEMENT DU NOMBRE MAXIMUM D'ENREGISTREMENT) LEUR
! --- PEUT ETRE NULLE
!
        call jxecro(ic, iadadi, iaddad, ladi, 0, &
                    lideff-1)
        iadd(jiadd(ic)+2*(lideff-1)-1) = iaddad(1)
        iadd(jiadd(ic)+2*(lideff-1)) = iaddad(2)
!
        iadad2 = iadm(jiadm(ic)+2*lideff-1)
        ladi2 = lono(jlono(ic)+lideff)*ltyp(jltyp(ic)+lideff)
        call jxecro(ic, iadad2, iaddac, ladi2, 0, &
                    lideff)
        iadd(jiadd(ic)+2*lideff-1) = iaddac(1)
        iadd(jiadd(ic)+2*lideff) = iaddac(2)
!
        idb = idebug
        idebug = 0
!
        do i = lideff-2, 2, -1
            iadmi = iadm(jiadm(ic)+2*i-1)
            if (iadmi .ne. 0) then
                idatos = i
                nomos = rnom(jrnom(ic)+i)
                call jjlide('SYSTEM', nomos, 1)
            end if
        end do
        do i = lidbas+1, nreuti(ic)
            iadmi = iadm(jiadm(ic)+2*i-1)
            iadyn = iadm(jiadm(ic)+2*i)
            call jjlidy(iadyn, iadmi)
        end do
        idebug = idb
    end if
!
    if (kcond .eq. 'SAUVE   ') then
!       ----------- STATISTIQUES DU FICHIER
!       ----------- ACTUALISER CARA
        cara(jcara(ic)+1) = nreuti(ic)
        cara(jcara(ic)+3) = nblmax(ic)
        cara(jcara(ic)+4) = nbluti(ic)
        cara(jcara(ic)+6) = iadd(jiadd(ic)+3)
        cara(jcara(ic)+7) = iadd(jiadd(ic)+4)
        cara(jcara(ic)+12) = nbenrg(ic)
        if (iadcar .ne. 0) then
            idatos = 1
            call jjlide('JELIBF', nomcar, 1)
            iadcar = 0
        end if
    end if
!
! ON PEUT MAINTENANT LIBERER LES ZONES MEMOIRES ASSOCIEES AUX
! OBJETS DU SYSTEME, ELLES NE SERONT PLUS UTILISEES (3-13):
!    $$GENR , $$TYPE , $$DOCU , $$ORIG , $$RNOM  , $$LTYP , $$LONG  ,
!    $$LONO , $$DATE , $$LUTI , $$HCOD
!
    if (kcond .ne. 'LIBERE  ') then
        do i = 3, lideff-2
            iadmi = iadm(jiadm(ic)+2*i-1)
            iadyn = iadm(jiadm(ic)+2*i)
            call jjlidy(iadyn, iadmi)
        end do
    end if
!
!
    if (kcond .eq. 'SAUVE   ') then
!       ----------- DECHARGER TAMPON D'ECRITURE ( PETITS OBJETS )
        lgbl = 1024*longbl(ic)*lois
        if (iitecr(ic) .gt. 0) then
            call jxecrb(ic, iitecr(ic), kitecr(ic)+1, lgbl, 0, &
                        0)
        end if
        if (litlec(ic)) then
            call jxecrb(ic, iitlec(ic), kitlec(ic)+1, lgbl, 0, &
                        0)
            litlec(ic) = .false.
            iitlec(ic) = 0
        end if
!       ---- ON DECHARGE MAINTENANT LES STATISTIQUES SUR LES ACCES
        idatos = lideff
        iacce(jiacce(ic)+iaddac(1)) = iacce(jiacce(ic)+iaddac(1))+1
        call jxecro(ic, iadacc, iaddac, lacc, 0, &
                    lideff)
        iitecr(ic) = 0
        if (litlec(ic)) then
            call jxecrb(ic, iitlec(ic), kitlec(ic)+1, lgbl, 0, &
                        0)
        end if
    end if
    if (kcond .eq. 'SAUVE   ' .or. kcond .eq. 'DETRUIT  ') then
!       ---- ON DECHARGE MAINTENANT LA DESCRIPTION DES ENREGISTREMENTS
        lgbl = 1024*longbl(ic)*lois
        idatos = lideff-1
        call jxecro(ic, iadadi, iaddad, ladi, 0, &
                    lideff-1)
        iitecr(ic) = 0
        if (litlec(ic)) then
            call jxecrb(ic, iitlec(ic), kitlec(ic)+1, lgbl, 0, &
                        0)
        end if
        call jjlidy(iadady, iadadi)
        iadady = 0
        iadadi = 0
!       ---- $$ACCE N'EST PLUS UTILISE, ON PEUT LE LIBERER
        call jjlidy(iadacy, iadacc)
        iadacy = 0
        iadacc = 0
    end if
!
    valk(1) = nombas(ic)
    vali(1) = nbluti(ic)
    vali(2) = nblmax(ic)
    vali(3) = 1024*longbl(ic)*lois
    vali(4) = nbacce(2*ic-1)
    valr(1) = nbacce(2*ic-1)*longbl(ic)*lois/1024.d0
    vali(5) = nbacce(2*ic)
    valr(2) = nbacce(2*ic)*longbl(ic)*lois/1024.d0
    vali(6) = nreuti(ic)
    vali(7) = nremax(ic)
    vali(8) = (nreuti(ic)*100)/nremax(ic)
    vali(9) = nbenrg(ic)
!
    if (info .ge. 1) then
        call utmess('I', 'JEVEUX_22', sk=valk(1), ni=9, vali=vali, &
                    nr=2, valr=valr)
    end if
!
    if (kcond .ne. 'LIBERE  ') then
        if (iadcar .ne. 0) then
            idatos = 1
            call jjlidy(iadcdy, iadcar)
        end if
        call jjlidy(iadady, iadacc)
!       ----------- CLORE LE FICHIER
        call jxferm(ic)
!       ----------- LIBERER PLACE
        call jjlidy(k18, k18i)
        call jjlidy(k2, iad2)
        call jjlidy(k19, k19i)
        call jjlidy(k16, kmarq(ic))
        call jjlidy(k17, k17i)
        call jjlidy(k20, kiadm(ic))
!
        if (kstout(ic) (1:7) .eq. 'DETRUIT') then
            nom = nomfic(ic) (1:4)//'.?  '
            call lxmins(nom)
            infrm = 0
            call rmfile(nom, infrm, iret)
        end if
!
        classe(ic:ic) = ' '
        nomfic(ic) = ' '
        nombas(ic) = ' '
        kstout(ic) = ' '
        kstini(ic) = ' '
        dn2(ic) = ' '
        nrhcod(ic) = 0
        nremax(ic) = 0
        nreuti(ic) = 0
        litlec(ic) = .false.
        nblmax(ic) = 0
        nbluti(ic) = 0
        longbl(ic) = 0
        kitlec(ic) = 0
        kitecr(ic) = 0
        kiadm(ic) = 0
        iitlec(ic) = 0
        iitecr(ic) = 0
        nitecr(ic) = 0
        jindir(ic) = 0
    end if
! FIN ------------------------------------------------------------------
end subroutine
