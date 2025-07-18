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
subroutine vefnme_cplx(option, base, &
                       model, mate, carael, &
                       comporZ, timePrev, timeCurr, nh, ligrelz, varicomz, &
                       sigmaPrev, sigmaz, strxz, deplz, vecelz)
!
    use HHO_precalc_module, only: hhoAddInputField
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dbgcal.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/codent.h"
#include "asterfort/exisd.h"
#include "asterfort/sepach.h"
#include "asterfort/copisd.h"
!
    character(len=16), intent(in) :: option
    character(len=1), intent(in) :: base
    character(len=8), intent(in) :: model
    real(kind=8), intent(in) :: timePrev, timeCurr
    character(len=8), intent(in) :: carael
    character(len=24), intent(in) :: mate
    character(len=*), intent(in) :: ligrelz
    integer(kind=8), intent(in) :: nh
    character(len=*), intent(in) :: comporZ
    character(len=*), intent(in) :: sigmaz, sigmaPrev
    character(len=*), intent(in) :: varicomz
    character(len=*), intent(in) :: strxz
    character(len=*), intent(in) :: deplz
    character(len=*), intent(inout) :: vecelz(*)
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Option: FORC_NODA
!         FONL_NOEU
!
! --------------------------------------------------------------------------------------------------
!
!
! IN  MODELE : NOM DU MODELE (NECESSAIRE SI SIGMA EST UNE CARTE)
! IN  SIGMA  : NOM DU CHAM_ELEM (OU DE LA CARTE) DE CONTRAINTES
! IN  CARA   : NOM DU CARA_ELEM
! IN  DEPMOI : NOM DU CHAM_NO DE DEPLACEMENTS PRECEDENTS
! IN  DEPDEL : NOM DU CHAM_NO D'INCREMENT DEPLACEMENTS
! IN  MATCOD : NOM DU MATERIAU CODE
! IN  COMPOR : NOM DE LA CARTE DE COMPORTEMENT
! IN  NH     : NUMERO D'HARMONIQUE DE FOURIER
! IN  PARTPS : INSTANT PRECEDENT ET ACTUEL
! IN  CARCRI : CARTE DES CRITERES ET DE THETA
! IN  CHVARC : NOM DU CHAMP DE VARIABLE DE COMMANDE
! IN  LIGREZ : (SOUS-)LIGREL DE MODELE POUR CALCUL REDUIT
!                  SI ' ', ON PREND LE LIGREL DU MODELE
! OUT VECELZ : VECT_ELEM RESULTAT.
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbout = 1
    integer(kind=8), parameter :: nbxin = 35
    character(len=8) :: lpaout(nbout), lpain(nbxin)
    character(len=19) :: lchout(nbout), lchin(nbxin)
!
    character(len=8) :: k8bla, mesh
    character(len=8) :: newnom, nomgd
    character(len=19) :: numhar, tpsmoi, tpsplu, ligrel_local, ligrel
    character(len=19) :: chgeom, chcara(18), vecele, veceli, compor
    character(len=19) :: lchinr(nbxin), lchini(nbxin)
    character(len=16) :: optio2
    integer(kind=8) :: iret, inddec(nbxin), iexi, k, nbin
    character(len=19) :: pintto, cnseto, heavto, loncha, basloc, lsn, lst, stano
    character(len=19) :: pmilto, fissno, hea_no
    character(len=19) :: sigma, varicom, strx
    character(len=19) :: depl
    character(len=19) :: chdecr(nbxin), chdeci(nbxin), ch19, chr, chi, ch1(nbout), ch2(nbout)
    aster_logical :: debug, lcmplx, lsspt
    integer(kind=8) :: ifmdbg, nivdbg
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! - Initializations
!
    sigma = sigmaz
    varicom = varicomz
    strx = strxz
    depl = deplz
    ligrel = ligrelz
    compor = comporZ
    newnom = '.0000000'
    numhar = '&&VEFNME.NUME_HARM'
    tpsmoi = '&&VEFNME.CH_INSTAM'
    tpsplu = '&&VEFNME.CH_INSTAP'
    k8bla = ' '
    optio2 = option
    if (option .ne. 'FONL_NOEU') optio2 = 'FORC_NODA'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
!
! - Get mesh
!
    if (depl .ne. ' ') then
        call dismoi('NOM_MAILLA', depl, 'CHAM_NO', repk=mesh)
    else if (sigma .ne. ' ') then
        call dismoi('NOM_MAILLA', sigma, 'CHAM_ELEM', repk=mesh)
    else
        ASSERT(.false.)
    end if
    chgeom = mesh(1:8)//'.COORDO'
!
! - VECT_ELEM name
!
    vecele = vecelz(1)
    if (vecele .eq. ' ') then
        vecele = '&&VEFNME'
    end if
    veceli = vecelz(2)
    if (veceli .eq. ' ') then
        veceli = '&&VEFNMI'
    end if
    if (ligrel .eq. ' ') then
        ligrel_local = model(1:8)//'.MODELE'
    else
        ligrel_local = ligrel
    end if
!
! - <CARTE> for structural elements
!
    call mecara(carael, chcara)
!
! - <CARTE> for Fourier mode
!
    call mecact('V', numhar, 'MAILLA', mesh, 'HARMON', &
                ncmp=1, nomcmp='NH', si=nh)
!
! - <CARTE> for instant
!
    call mecact('V', tpsmoi, 'MAILLA', mesh, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=timePrev)
    call mecact('V', tpsplu, 'MAILLA', mesh, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=timeCurr)
!
! - Init fields
    call inical(nbxin, lpain, lchin, nbout, lpaout, lchout)
    call inical(nbxin, lpain, lchin, nbout, lpaout, ch1)
    call inical(nbxin, lpain, lchin, nbout, lpaout, ch2)
!
! - CREATION DES LISTES DES CHAMPS IN
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = mate(1:19)
    lpain(3) = 'PCAGNPO'
    lchin(3) = chcara(6)
    lpain(4) = 'PCAORIE'
    lchin(4) = chcara(1)
    lpain(5) = 'PCOMPOR'
    lchin(5) = compor
    lpain(6) = 'PSIEFR'
    lchin(6) = sigmaz
    lpain(7) = 'PDEPLAR'
    lchin(7) = depl
    lpain(8) = 'PCONTGM'
    lchin(8) = sigmaPrev
    lpain(9) = 'PCAARPO'
    lchin(9) = chcara(9)
    lpain(10) = 'PCADISK'
    lchin(10) = chcara(2)
    lpain(11) = 'PCACOQU'
    lchin(11) = chcara(7)
    lpain(12) = 'PHARMON'
    lchin(12) = numhar
    lpain(13) = 'PCAMASS'
    lchin(13) = chcara(12)
    lpain(14) = 'PINSTMR'
    lchin(14) = tpsmoi
    lpain(15) = 'PINSTPR'
    lchin(15) = tpsplu
    lpain(16) = 'PVARCPR'
    lchin(16) = varicom
    lpain(17) = 'PCAGEPO'
    lchin(17) = chcara(5)
    lpain(18) = 'PNBSP_I'
    lchin(18) = chcara(16)
    lpain(19) = 'PFIBRES'
    lchin(19) = chcara(17)

! --- CADRE X-FEM
    call exixfe(model, iret)
    if (iret .ne. 0) then
        pintto = model(1:8)//'.TOPOSE.PIN'
        cnseto = model(1:8)//'.TOPOSE.CNS'
        heavto = model(1:8)//'.TOPOSE.HEA'
        loncha = model(1:8)//'.TOPOSE.LON'
        pmilto = model(1:8)//'.TOPOSE.PMI'
        basloc = model(1:8)//'.BASLOC'
        lsn = model(1:8)//'.LNNO'
        lst = model(1:8)//'.LTNO'
        stano = model(1:8)//'.STNO'
        fissno = model(1:8)//'.FISSNO'
        hea_no = model(1:8)//'.TOPONO.HNO'
    else
        pintto = '&&VEFNME.PINTTO.BID'
        cnseto = '&&VEFNME.CNSETO.BID'
        heavto = '&&VEFNME.HEAVTO.BID'
        loncha = '&&VEFNME.LONCHA.BID'
        basloc = '&&VEFNME.BASLOC.BID'
        pmilto = '&&VEFNME.PMILTO.BID'
        lsn = '&&VEFNME.LNNO.BID'
        lst = '&&VEFNME.LTNO.BID'
        stano = '&&VEFNME.STNO.BID'
        fissno = '&&VEFNME.FISSNO.BID'
        hea_no = '&&VEFNME.HEA_NO.BID'
    end if
!
    lpain(20) = 'PPINTTO'
    lchin(20) = pintto
    lpain(21) = 'PCNSETO'
    lchin(21) = cnseto
    lpain(22) = 'PHEAVTO'
    lchin(22) = heavto
    lpain(23) = 'PLONCHA'
    lchin(23) = loncha
    lpain(24) = 'PBASLOR'
    lchin(24) = basloc
    lpain(25) = 'PLSN'
    lchin(25) = lsn
    lpain(26) = 'PLST'
    lchin(26) = lst
    lpain(27) = 'PSTANO'
    lchin(27) = stano
    lpain(28) = 'PCINFDI'
    lchin(28) = chcara(15)
    lpain(29) = 'PPMILTO'
    lchin(29) = pmilto
    lpain(30) = 'PFISNO'
    lchin(30) = fissno
    lpain(31) = 'PSTRXMR'
    lchin(31) = strx
    lpain(32) = 'PHEA_NO'
    lchin(32) = hea_no
!
    nbin = 32
!
    call hhoAddInputField(model, nbxin, lchin, lpain, nbin)
!
! --- CREATION DES LISTES DES CHAMPS OUT
!
    lpaout(1) = 'PVECTUR'
    call gcnco2(newnom)
    lchout(1) = vecele(1:8)//newnom
    call corich('E', lchout(1), ichin_=-1)
    call gcnco2(newnom)
    ch1(1) = vecele(1:8)//newnom
    call corich('E', ch1(1), ichin_=-1)
    call gcnco2(newnom)
    ch2(1) = veceli(1:8)//newnom
    call corich('E', ch2(1), ichin_=-1)
!
! --- PREPARATION DU VECT_ELEM RESULTAT
!
    call detrsd('VECT_ELEM', vecele)
    call memare(base, vecele, model, 'CHAR_MECA')
!
    lcmplx = .false.
    do k = 1, nbin
        inddec(k) = 0
        if (lpain(k) .eq. ' ') goto 1
        ch19 = lchin(k)
        if (ch19 .eq. ' ') goto 1
        call exisd('CHAMP', ch19, iexi)
        if (iexi .eq. 0) goto 1
        call dismoi('NOM_GD', ch19, 'CHAMP', repk=nomgd)
        if (nomgd(5:6) .eq. '_C') then
            lcmplx = .true.
            inddec(k) = 1
            chr = '&&VEFNME.CHXX.R'
            chi = '&&VEFNME.CHXX.I'
            call codent(k, 'D0', chr(12:13))
            call codent(k, 'D0', chi(12:13))
            call sepach(carael, ch19, 'V', chr, chi)
            chdecr(k) = chr
            chdeci(k) = chi
        end if
1       continue
    end do

    if (lcmplx) then
        call detrsd('VECT_ELEM', veceli)
    end if

    call exisd('CHAM_ELEM_S', lchout(1), iexi)
    lsspt = (iexi .ne. 0)

    if (lcmplx) then
        do k = 1, nbin
            if (inddec(k) .eq. 0) then
                lchinr(k) = lchin(k)
                lchini(k) = lchin(k)
            else
                lchinr(k) = chdecr(k)
                lchini(k) = chdeci(k)
            end if
        end do
        call memare(base, veceli, model, 'CHAR_MECA')
    end if

    if (debug) then
        call dbgcal(optio2, ifmdbg, nbin, lpain, lchin, &
                    nbout, lpaout, lchout)
    end if
!
! - APPEL A CALCUL
!
    if (lcmplx) then
        if (lsspt) call copisd('CHAM_ELEM_S', 'V', lchout(1), ch1(1))
        call calcul('S', optio2, ligrel_local, nbin, lchinr, &
                    lpain, nbout, ch1, lpaout, 'V', &
                    'OUI')
        call reajre(vecele, ch1(1), 'V')
        vecelz(1) = vecele//'.RELR'
!
        if (lsspt) call copisd('CHAM_ELEM_S', 'V', lchout(1), ch2(1))
        call calcul('S', optio2, ligrel_local, nbin, lchini, &
                    lpain, nbout, ch2, lpaout, 'V', &
                    'OUI')
        call reajre(veceli, ch2(1), 'V')
        vecelz(2) = veceli//'.RELR'
    else
        call calcul('S', optio2, ligrel_local, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
        call reajre(vecele, lchout(1), base)
        vecelz(1) = vecele//'.RELR'
    end if
!
! - Cleaning
!
    do k = 1, nbin
        if (inddec(k) .ne. 0) then
            call detrsd('CHAMP', chdecr(k))
            call detrsd('CHAMP', chdeci(k))
        end if
    end do

    call detrsd('CHAMP_GD', numhar)
    call jedema()
end subroutine
