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
subroutine dmgmod(nomsym, nomsd, nomsd2, nommat, nbordr, &
                  jordr, jcoef, nbpt, ntcmp, numcmp, &
                  impr, vdomag)
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rccome.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nommat
    character(len=16) :: nomsym
    character(len=19) :: nomsd, nomsd2
    real(kind=8) :: vdomag(*)
    integer(kind=8) :: nbpt, numcmp(*)
    integer(kind=8) :: ntcmp, impr, nbordr, jordr, jcoef
!       CREATION D UN VECTEUR AUX NOEUDS/PG : AMPLITUDE MAX DE VIBRATION
!       METHODE CALCUL DU DOMMAGE UNITAIRE = /WOHLER
!       ----------------------------------------------------------------
!       IN     NOMSYM    NOM SYMBOLIQUE OPTION EQUI_GD
!              NOMSD     NOM SD RESULTAT STATIQUE
!              NOMSD2     NOM SD RESULTAT MODAL
!              NOMMAT    NOM DU CHAM_MATER
!              NBORDR    NOMBRE DE NUMEROS D'ORDRE
!              JORD      ADRESSE DE LA LISTE DES NUMEROS D'ORDRE
!              JCEIF      ADRESSE DE LA LISTE DES COEFFICIENTS
!              NBPT      NOMBRE DE POINTS DE CALCUL DU DOMMAGE
!              NTCMP     NOMBRE TOTAL DE COMPOSANTE OPTION EQUI_GD
!              NUMCMP    NUMERO(S) DE LA(DES) COMPOSANTE(S) DE EQUI_GD
!              IMPR      NIVEAU IMPRESSION
!       OUT    VDOMAG    VECTEUR DOMMAGE AUX POINTS
!       ----------------------------------------------------------------
!       ---------------------------------------------------------------
    character(len=11) :: k11
    character(len=8) :: nompar, kcorre, nomfon
    character(len=16) :: nomrm(1)
    character(len=32) :: nomphe
    character(len=19) :: chequi, chequ2(nbordr)
    character(len=24) :: nomdmg
    character(len=24) :: valk(3)
    integer(kind=8) :: icodre(1)
!
    real(kind=8) :: su, salt0, dmax, saltm, val(1)
    real(kind=8) :: valr(3), r8b, dmin, smax, coeff, r8min
!
    integer(kind=8) :: ipt, iord, icmp, nbr, nbk, nbc, nbf
    integer(kind=8) :: ivch, ivpt, ibid, ivalk
    integer(kind=8) :: numsym, ivch2, ivord2(nbordr), numord
    aster_logical :: crit
!
! ---   VECTEURS DE TRAVAIL
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ik
    real(kind=8), pointer :: celv(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    r8b = 0.d0
    call jemarq()
!
    nomdmg = '&&OP0151.EQUI_GD'
    call wkvect(nomdmg, 'V V R', 2, ivpt)
!
    call getvtx(' ', 'CORR_SIGM_MOYE', scal=kcorre, nbret=ibid)
!
! --    VECTEUR DES NUMORD NOMS DE CHAMPS POUR L OPTION NOMSYM
!
    call jenonu(jexnom(nomsd//'.DESC', nomsym), numsym)
    if (numsym .eq. 0) then
        valk(1) = nomsym
        valk(2) = nomsd
        call utmess('F', 'PREPOST_51', nk=2, valk=valk)
    end if
    call jeveuo(jexnum(nomsd//'.TACH', numsym), 'L', ivch)
!
    call jenonu(jexnom(nomsd2//'.DESC', nomsym), numsym)
    if (numsym .eq. 0) then
        valk(1) = nomsym
        valk(2) = nomsd2
        call utmess('F', 'PREPOST_51', nk=2, valk=valk)
    end if
    call jeveuo(jexnum(nomsd2//'.TACH', numsym), 'L', ivch2)
!
! RECUPERATION PROPRIETES MATERIAUX
!
    nomrm(1) = 'SU'
    nompar = ' '
    call rcvale(nommat, 'RCCM', 0, nompar, [r8b], &
                1, nomrm, val, icodre(1), 2)
    if (icodre(1) .ne. 0) then
        valk(1) = 'SU'
        call utmess('F', 'FATIGUE1_88', sk=valk(1))
    end if
    su = val(1)
!
    nomphe = 'FATIGUE   '
    call rccome(nommat, nomphe, icodre(1), k11_ind_nomrc=k11)
    call jelira(nommat//k11//'.VALR', 'LONUTI', nbr)
    call jelira(nommat//k11//'.VALC', 'LONUTI', nbc)
    call jeveuo(nommat//k11//'.VALK', 'L', ivalk)
    call jelira(nommat//k11//'.VALK', 'LONUTI', nbk)
    nbf = (nbk-nbr-nbc)/2
    do ik = 1, nbf
        if (zk16(ivalk-1+nbr+nbc+ik) .eq. 'WOHLER') then
            nomfon = zk16(ivalk-1+nbr+nbc+nbf+ik)
            call jeveuo(nomfon//'           .VALE', 'L', vr=vale)
            salt0 = vale(1)
        end if
    end do
!
    valr(1) = su
    valr(2) = salt0
    call utmess('I', 'FATIGUE1_87', nr=2, valr=valr)
!
    icmp = 1
    dmin = 1.d10
    smax = 0.d0
    crit = .false.
    r8min = r8miem()
!
! ---       CALCUL DU VECTEUR HISTOIRE DE LA EQUI_GD EN CE POINT
!
    chequi = zk24(ivch) (1:19)
    if (chequi .eq. ' ') then
        valk(1) = chequi
        valk(2) = nomsym
        valk(3) = nomsd
        call utmess('F', 'PREPOST_52', nk=3, valk=valk)
    end if
    call jeveuo(chequi//'.CELV', 'L', vr=celv)
!
    do iord = 1, nbordr
        numord = zi(jordr+iord-1)
        chequ2(iord) = zk24(ivch2+numord-1) (1:19)
        if (chequ2(iord) .eq. ' ') then
            valk(1) = chequ2(iord)
            valk(2) = nomsym
            valk(3) = nomsd2
            call utmess('F', 'PREPOST_52', nk=3, valk=valk)
        end if
        call jeveuo(chequ2(iord)//'.CELV', 'L', ivord2(iord))
    end do
!
! ---     BOUCLE SUR LES POINTS
!
    do ipt = 1, nbpt
! -    STOCKAGE CONTRAINTES
        zr(ivpt) = celv(1+(ipt-1)*ntcmp+numcmp(icmp)-1)
        zr(ivpt+1) = 0.d0
        do iord = 1, nbordr
            coeff = zr(jcoef+iord-1)
            zr(ivpt+1) = zr(ivpt+1)+coeff*abs(zr(ivord2(iord)+(ipt-1)*ntcmp+numcmp(icmp)-1))
        end do
!
        if (zr(ivpt) .gt. su) then
            if (impr .ge. 2) then
                call utmess('I', 'FATIGUE1_80')
            end if
            saltm = 0.d0
            crit = .true.
            smax = max(zr(ivpt), smax)
        else if (zr(ivpt) .gt. 0.d0) then
            if (kcorre .eq. 'GOODMAN') then
                saltm = salt0*(1-zr(ivpt)/su)
            else if (kcorre .eq. 'GERBER') then
                saltm = salt0*(1-(zr(ivpt)/su)**2)
            end if
        else
            saltm = salt0
        end if
!
        if (abs(zr(ivpt+1)) .gt. r8min) then
            dmax = saltm/abs(zr(ivpt+1))
        else
            dmax = 1.d10
        end if
        if (dmax .lt. dmin) dmin = dmax
!
        vdomag(ipt) = dmax
        if (impr .ge. 2) then
            valr(1) = zr(ivpt)
            valr(2) = zr(ivpt+1)
            valr(3) = dmax
            call utmess('I', 'FATIGUE1_79', si=ipt, nr=3, valr=valr)
        end if
!
    end do
!
    if (crit) then
        valr(1) = smax
        valr(2) = su
        call utmess('A', 'FATIGUE1_83', nr=2, valr=valr)
    end if
    call utmess('I', 'FATIGUE1_82', sr=dmin)
!
!
    call jedema()
end subroutine
