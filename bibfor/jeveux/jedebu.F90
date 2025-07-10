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

subroutine jedebu(nbfi, mxzon, idb)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "jeveux_private.h"
#include "asterc/gtoptk.h"
#include "asterc/gtoptr.h"
#include "asterc/ismaem.h"
#include "asterc/isnnem.h"
#include "asterc/ispbem.h"
#include "asterc/lbisem.h"
#include "asterc/loc8em.h"
#include "asterc/lofiem.h"
#include "asterc/loisem.h"
#include "asterc/lolsem.h"
#include "asterc/lor8em.h"
#include "asterc/mofiem.h"
#include "asterc/r8nnem.h"
#include "asterfort/adjust_memlimit.h"
#include "asterfort/assert.h"
#include "asterfort/jxdate.h"
#include "asterfort/utmess.h"
#include "asterfort/utptme.h"
    integer(kind=8) :: nbfi, mxzon, idb
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR D'INITIALISATION GENERALE POUR LE GESTIONNAIRE
!         DE MEMOIRE
!
! IN  NBFI   : NOMBRE MAXIMUM DE BASES SIMULTANEES ( =< 5)
! IN  MXZON  : LIMITE MEMOIRE DYNAMIQUE (EN ENTIER words)
! IN  IDB    : PARAMETRE DE DEBUG
!              ( 0: RIEN, 1: MISE A UNDEF DES SEGMENTS DE VALEURS )
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    integer(kind=8) :: nbfic
    common/iparje/nbfic
    integer(kind=8) :: iloc
    common/ilocje/iloc
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, jcara, jdate, jdocu, jgenr, jhcod
    integer(kind=8) :: jiacce, jiadd, jiadm, jlong, jlono, jltyp
    integer(kind=8) :: jluti, jmarq, jorig, jrnom, jtype, k
    integer(kind=8) :: n, nbacce
    real(kind=8) :: val
!-----------------------------------------------------------------------
    parameter(n=5)
!
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    aster_logical :: litlec
    common/lficje/litlec(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
! ----------------------------------------------------------------------
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
! ----------------------------------------------------------------------
    integer(kind=8) :: isstat
    common/iconje/isstat
    integer(kind=8) :: msstat, lsstat
    common/jconje/msstat, lsstat
! ----------------------------------------------------------------------
    integer(kind=8) :: datei
    common/iheuje/datei
! ----------------------------------------------------------------------
    integer(kind=8) :: illici, jclass(0:255)
    common/jchaje/illici, jclass
! ----------------------------------------------------------------------
    integer(kind=8) :: istat
    common/istaje/istat(4)
    character(len=4) :: kstat
    common/kstaje/kstat
    integer(kind=8) :: mslois
    common/jenvje/mslois
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: idn, iext, nbenrg
    common/iextje/idn(n), iext(n), nbenrg(n)
    integer(kind=8) :: lfic, mfic
    common/fenvje/lfic(n), mfic
    character(len=128) :: repglo, repvol
    common/banvje/repglo, repvol
    integer(kind=8) :: lrepgl, lrepvo
    common/balvje/lrepgl, lrepvo
    integer(kind=8) :: lundef, idebug
    common/undfje/lundef, idebug
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    integer(kind=8) :: icdyn, mxltot
    common/xdynje/icdyn, mxltot
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
    common/jiacce/jiacce(n), nbacce(2*n)
    integer(kind=8) :: indiq_jjagod, indiq_jjldyn
    common/idagod/indiq_jjagod, indiq_jjldyn

! --------------------------------- ------------------------------------
    integer(kind=8) :: mxlici, iret
    parameter(mxlici=67)
    character(len=mxlici) :: clicit
    data clicit/' ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.$&_abcdefghijklmnopqrstuvwxyz'/
! DEB ------------------------------------------------------------------
!
! ON AFFECTE ZI(1), ZR(1) ET ZC(1) AVEC UNE VALEUR QUI PEUT FAIRE
! PLANTER SI ELLE EST UTILISEE
!
    zi(1) = ismaem()
    zr(1) = r8nnem()
    zc(1) = dcmplx(zr(1), zr(1))
! -----------------  ENVIRONNEMENT MACHINE -----------------------------
    indiq_jjldyn = 0
    indiq_jjagod = 0
    do k = 1, n
        lfic(k) = lofiem()
    end do
    call gtoptr('maxbase', val, iret)
    if (val .le. 0 .or. iret .ne. 0) then
        mfic = mofiem()
    else
        mfic = nint(val)*1024
    end if
    call gtoptk('repglob', repglo, iret)
    if (iret .ne. 0) then
        repglo = '. '
        lrepgl = 1
    else
        lrepgl = index(repglo, ' ')-1
        if (lrepgl .gt. 119) then
            call utmess('F', 'JEVEUX1_69', sk=repglo, si=lrepgl)
        end if
    end if
    call gtoptk('repvola', repvol, iret)
    if (iret .ne. 0) then
        repvol = '. '
        lrepvo = 1
    else
        lrepvo = index(repvol, ' ')-1
        if (lrepvo .gt. 119) then
            call utmess('F', 'JEVEUX1_70', sk=repvol, si=lrepvo)
        end if
    end if
    lbis = lbisem()
    lor8 = lor8em()
    loc8 = loc8em()
    lois = loisem()
    lols = lolsem()
    lundef = isnnem()
    mslois = lois-1
    ldyn = 1
    lgdyn = 1
    mxdyn = 0
    lgio(1) = 0
    lgio(2) = 0
    mcdyn = 0
    mldyn = 0
    nbdyn = 0
    nbfree = 0
    icdyn = 0
    mxltot = 0
    svuse = 16
    smxuse = svuse
! -----------------  NOMBRE DE BASES -----------------------------------
    nbfic = min(nbfi, n, len(classe))
    ASSERT(nbfic .gt. 0 .and. nbfic .eq. nbfi)
! -----------------  CONSTANTES DE STATUT DES SEGMENTS DE VALEURS ------
    kstat = 'XUAD'
    isstat = ispbem(lbis-3)
    do k = 1, 4
        istat(k) = k*isstat
    end do
    idebug = idb
! -----------------  ZONE MEMOIRE  -------------------------------------
    vmxdyn = mxzon
    if (mxzon .eq. 0) then
        vmxdyn = 1024
    end if
!   memoire totale demandee par l'utilisateur
    vmet = vmxdyn
!
    call utptme('MEM_MUMP', 0.d0, iret)
!   ajustement pour démarrer, un second sera fait dans debut après lecture des catalogues
    call adjust_memlimit(ASTER_FALSE)

    liszon = 1
    jiszon = 1
    lk1zon = liszon*lois
    jk1zon = jiszon*lois
    iloc = loc(iszon(jiszon))
! -------------------  POINTEURS D'ATTRIBUTS  --------------------------
    do i = 1, len(classe)
        classe(i:i) = '$'
    end do
    do i = 1, nbfic
        jgenr(i) = 0
        jtype(i) = 0
        jltyp(i) = 0
        jdocu(i) = 0
        jorig(i) = 0
        jrnom(i) = 0
        jlono(i) = 0
        jlong(i) = 0
        jdate(i) = 0
        jiadd(i) = 0
        jiadm(i) = 0
        jmarq(i) = 0
        jluti(i) = 0
        jcara(i) = 0
        jhcod(i) = 0
        nremax(i) = 0
        nrhcod(i) = 0
        nreuti(i) = 0
        nblmax(i) = 0
        nbenrg(i) = 1
        nbluti(i) = 0
        longbl(i) = 0
        kitlec(i) = 0
        kitecr(i) = 0
        kiadm(i) = 0
        iitlec(i) = 0
        iitecr(i) = 0
        nitecr(i) = 0
        litlec(i) = .false.
        nomfic(i) = '        '
        kstout(i) = '        '
        kstini(i) = '        '
        classe(i:i) = ' '
        dn2(i) = ' '
        nbacce(2*i-1) = 0
        nbacce(2*i) = 0
    end do
! -------------------  CONSTANTES DE GESTION  --------------------------
    lsstat = lbis-4
    msstat = 0
    do k = 1, lbis-4
        msstat = msstat+ispbem(k)
    end do
    bl32 = ' '
!
    call jxdate(datei)
!
    illici = -1
    do k = 0, 255
        jclass(k) = illici
    end do
    do k = 1, mxlici
        jclass(ichar(clicit(k:k))) = k
    end do
!
    kdesma(1) = 0
    kdesma(2) = 0
    lgduti = 0
    kposma(1) = 0
    kposma(2) = 0
    lgputi = 0
    ipgc = 0
! FIN ------------------------------------------------------------------
end subroutine
