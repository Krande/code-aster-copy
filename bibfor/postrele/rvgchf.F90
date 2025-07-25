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
subroutine rvgchf(epsi, criter, nomsd, chpsym, acces, &
                  ival, rval, nbval, ncheff)
    implicit none
!
#include "jeveux.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: criter
    character(len=16) :: chpsym, ncheff
    character(len=8) :: nomsd
    character(len=2) :: acces
    integer(kind=8) :: ival(*), nbval
    real(kind=8) :: rval(*), epsi
!
!     GENERATION DE LA LISTE DES NOM DE CHAMP EFFECTIFS DANS UNE
!     SD RESULTAT POUR UN CHAMP SYMBOLIQUE
!     ------------------------------------------------------------------
! IN  EPSI   : R : PRECISION DEMANDEE
! IN  CRITER : K : CRITERE DE COMPARAISON DE DEUX REELS
! IN  NOMSD  : K : NOM DE LA SD RESULTAT
! IN  CHPSYM : K : NOM DU CHAMP SYMBOLIQUE
! IN  ACCES  : K : TYPE D' ACCES DEMANDE
! IN  IVAL   : I : TABLE DES VALEURS ENTIERES POUR L' ACCES
! IN  RVAL   : R : TABLE DES VALEURS REELLES POUR L' ACCES
! IN  NBVAL  : I : DIMENSION DE TABLES XVAL
! OUT NCHEFF : K : NOM DE L' OJB DES NOM DE CHAMPS EFFECTIFS
!     ------------------------------------------------------------------
!     NCHEFF ::= RECORD
!       '.TYPACCE' : V E K8
!                    'INSTANT ', 'FREQUENC', 'MODE    ', 'ORDRE   '
!                    OU 'DIRECT   ' QUI CORRESPOND AU CAS PARTICULIER
!                    DU CHAMP_GD (NON TRAITE ICI MAIS DANS OP0051)
!       '.VALACCE' : V V SCAL
!                    CONTIENT LES VALEURS UTILISEES POUR L' ACCES
!       '.LSCHEFF' : XD V V K24 NUMEROTEE
!                    LES NOM DE CHAMPS EFFECTIFS CORRESPONDANT A LA
!                    VALEUR NUMERO I SONT RANGES DANS L' OC NUMERO I
!                    '&...' CODE LA NON PRESENCE D' UN CHAMP EFFECTIF
!     ------------------------------------------------------------------
!
!
!
    character(len=24) :: ntypac, nvalac, nlschp, nnores
    character(len=24) :: valk
    character(len=16) :: modacc
    character(len=8) :: k8bid
    character(len=1) :: car1, car2
    integer(kind=8) :: n1, n2, n3, i, j, nbordr, ibid, iordr, nbtrou
    integer(kind=8) :: aliste, atypac, avalac, alschp, avalr8, avalis, anores
    integer(kind=8) :: vali, tord(1)
    real(kind=8) :: rbid
    complex(kind=8) :: cbid
!
!======================================================================
!
    call jemarq()
    car1 = acces(1:1)
    car2 = acces(2:2)
!
    nnores = ncheff//'.NOMRESU'
    ntypac = ncheff//'.TYPACCE'
    nvalac = ncheff//'.VALACCE'
    nlschp = ncheff//'.LSCHEFF'
!
    call wkvect(ntypac, 'V V K8', 1, atypac)
    call wkvect(nnores, 'V V K16', 2, anores)
    zk16(anores) = nomsd
    zk16(anores+1) = chpsym
!
    if (car1 .eq. 'T') then
!     /* ACCES A TOUS LES NUMEROS D' ORDRES */
        cbid = dcmplx(0, 0)
        call rsorac(nomsd, 'LONUTI', ibid, rbid, k8bid, &
                    cbid, rbid, k8bid, tord, 1, &
                    nbtrou)
        nbordr = tord(1)
        call wkvect(nvalac, 'V V I', nbordr, avalac)
        call rsorac(nomsd, 'TOUT_ORDRE', ibid, rbid, k8bid, &
                    cbid, rbid, k8bid, zi(avalac), nbordr, &
                    nbtrou)
        zk8(atypac) = 'ORDRE   '
        call jecrec(nlschp, 'V V K24', 'NU', 'DISPERSE', 'VARIABLE', &
                    nbordr)
        do i = 1, nbordr, 1
            iordr = zi(avalac+i-1)
            call jecroc(jexnum(nlschp, i))
            call jeecra(jexnum(nlschp, i), 'LONMAX', 1)
            call jeveuo(jexnum(nlschp, i), 'E', alschp)
            call rsexch(' ', nomsd, chpsym, iordr, zk24(alschp+1-1), &
                        n1)
            if (n1 .ne. 0) then
                valk = chpsym
                vali = iordr
                call utmess('A', 'POSTRELE_41', sk=valk, si=vali)
                zk24(alschp+1-1) = '&&CHAMP_EFF_NON_EXISTANT'
            end if
        end do
    else
!     /* ACCES PAR LISTES ENUMEREES */
        if (car1 .ne. 'I') then
            nbordr = nbval
            if ((car2 .eq. 'O') .or. (car2 .eq. 'M')) then
                call wkvect('&&OP0051.LISTE.IS', 'V V I', nbordr, avalis)
                do i = 1, nbordr, 1
                    zi(avalis+i-1) = ival(i)
                end do
            else
                call wkvect('&&OP0051.LISTE.R8', 'V V R', nbordr, avalr8)
                do i = 1, nbordr, 1
                    zr(avalr8+i-1) = rval(i)
                end do
            end if
        end if
        call jecrec(nlschp, 'V V K24', 'NU', 'DISPERSE', 'VARIABLE', &
                    nbordr)
        if (car2 .eq. 'O') then
!        /* CAS D' UNE LISTE DE NUMERO ORDRE */
            call wkvect(nvalac, 'V V I', nbordr, avalac)
            zk8(atypac) = 'ORDRE   '
            do i = 1, nbordr, 1
                zi(avalac+i-1) = zi(avalis+i-1)
            end do
            do j = 1, nbordr, 1
                call jecroc(jexnum(nlschp, j))
                call jeecra(jexnum(nlschp, j), 'LONMAX', 1)
                call jeveuo(jexnum(nlschp, j), 'E', alschp)
                call rsexch(' ', nomsd, chpsym, zi(avalac+j-1), zk24(alschp+1-1), &
                            n2)
                if (n2 .ne. 0) then
                    valk = chpsym
                    vali = zi(avalac+j-1)
                    call utmess('A', 'POSTRELE_41', sk=valk, si=vali)
                    zk24(alschp+1-1) = '&&CHAMP_EFF_NON_EXISTANT'
                end if
            end do
        else if (car2 .eq. 'M') then
!        /* CAS D' UNE LISTE DE NUMERO DE MODE */
            call wkvect(nvalac, 'V V I', nbordr, avalac)
            zk8(atypac) = 'MODE    '
            modacc = 'NUME_MODE'
            do i = 1, nbordr, 1
                zi(avalac+i-1) = zi(avalis+i-1)
            end do
            do i = 1, nbordr, 1
                cbid = dcmplx(0, 0)
                call rsorac(nomsd, modacc, zi(avalac+i-1), 0.d0, k8bid, &
                            cbid, epsi, criter, zi, 0, &
                            n1)
                n1 = -n1
                call jecroc(jexnum(nlschp, i))
                n3 = max(n1, 1)
                call jeecra(jexnum(nlschp, i), 'LONMAX', n3)
                call jeveuo(jexnum(nlschp, i), 'E', alschp)
                if (n1 .eq. 0) then
                    zk24(alschp+1-1) = '&&CHAMP_EFF_NON_EXISTANT'
                else
                    call wkvect('&&OP0051.LISTE.ORDRE', 'V V I', n1, aliste)
                    cbid = dcmplx(0, 0)
                    call rsorac(nomsd, modacc, zi(avalac+i-1), 0.0d0, k8bid, &
                                cbid, epsi, criter, zi(aliste), n1, &
                                n2)
                    do j = 1, n1, 1
                        call rsexch(' ', nomsd, chpsym, zi(aliste+j-1), zk24(alschp+j-1), &
                                    n2)
                        if (n2 .ne. 0) then
                            valk = chpsym
                            vali = zi(aliste+j-1)
                            call utmess('A', 'POSTRELE_41', sk=valk, si=vali)
                            zk24(alschp+j-1) = '&&CHAMP_EFF_NON_EXISTANT'
                        end if
                    end do
                    call jedetr('&&OP0051.LISTE.ORDRE')
                end if
            end do
        else
!        /* CAS D' UNE LISTE DE REELS */
            call wkvect(nvalac, 'V V R8', nbordr, avalac)
            if (car2 .eq. 'I') then
                zk8(atypac) = 'INSTANT '
                modacc = 'INST'
            else
                zk8(atypac) = 'FREQUENC'
                modacc = 'FREQ'
            end if
            do i = 1, nbordr, 1
                zr(avalac+i-1) = zr(avalr8+i-1)
            end do
            do i = 1, nbordr, 1
                cbid = dcmplx(0, 0)
                call rsorac(nomsd, modacc, 0, zr(avalac+i-1), k8bid, &
                            cbid, epsi, criter, zi, 0, &
                            n1)
                n1 = -n1
                call jecroc(jexnum(nlschp, i))
                n3 = max(n1, 1)
                call jeecra(jexnum(nlschp, i), 'LONMAX', n3)
                call jeveuo(jexnum(nlschp, i), 'E', alschp)
                if (n1 .eq. 0) then
                    zk24(alschp+1-1) = '&&CHAMP_EFF_NON_EXISTANT'
                    call utmess('A', 'POSTRELE_49', sr=zr(avalac+i-1))
                else
                    call jecreo('&&OP0051.LISTE.ORDRE', 'V V I')
                    call jeecra('&&OP0051.LISTE.ORDRE', 'LONMAX', n1)
                    call jeveuo('&&OP0051.LISTE.ORDRE', 'E', aliste)
                    cbid = dcmplx(0, 0)
                    call rsorac(nomsd, modacc, 0, zr(avalac+i-1), k8bid, &
                                cbid, epsi, criter, zi(aliste), n1, &
                                n2)
                    do j = 1, n1, 1
                        call rsexch(' ', nomsd, chpsym, zi(aliste+j-1), zk24(alschp+j-1), &
                                    n2)
                        if (n2 .ne. 0) then
                            valk = chpsym
                            vali = zi(aliste+j-1)
                            call utmess('A', 'POSTRELE_41', sk=valk, si=vali)
                            zk24(alschp+j-1) = '&&CHAMP_EFF_NON_EXISTANT'
                        end if
                    end do
                    call jedetr('&&OP0051.LISTE.ORDRE')
                end if
            end do
        end if
        call jedetr('&&OP0051.LISTE.IS')
        call jedetr('&&OP0051.LISTE.R8')
    end if
    call jedema()
end subroutine
