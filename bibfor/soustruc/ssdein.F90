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
subroutine ssdein(chno_lz, chno_gz, mail, nocas)
    implicit none
!     ARGUMENTS:
!     ----------
! ----------------------------------------------------------------------
!     BUT:
!      - CALCULER LE CHAMP DE DEPLACEMENT INTERNE A UNE SOUS-STRUCTURE
!        A PARTIR DU CHAMP DE DEPLACEMENT CONNU SUR SES NOEUDS EXTERNES
!
! IN_F,OU_J: CHNO_L : NOM DU CHAMP LOCAL A LA SOUS-STRUCTURE
! IN_F,IN_J: CHNO_G : NOM DU CHAMP GLOBAL (MODELE DE NIVEAU SUPERIEUR)
! IN_F     : MAIL : NOM DE LA (SUPER)MAILLE SUR LAQUELLE ON VEUT CHNO_L
! IN_F     : NOCAS: NOM DU CHARGEMENT CORRESPONDANT (EN PRINCIPE) A CHNO_G.
!                   (EVENTUELLEMENT : ' ')
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/matrot.h"
#include "asterfort/nueq_chck.h"
#include "asterfort/ssrone.h"
#include "asterfort/ssvaro.h"
#include "asterfort/ssvau1.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: chno_lz, chno_gz, mail, nocas, mag, mal, nomgd, nomacr
    character(len=14) :: nul
    character(len=19) :: nume_equa_g, nume_equa_l
    real(kind=8) :: lambda(6, 6), angl(3), pgl(3, 3)
    aster_logical :: exil, exig
    character(len=8) :: rota, ch8(2)
    character(len=19) :: chno_g, chno_l
    character(len=24) :: valk(2)
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadgg, iadgl, ialica
    integer(kind=8) :: ialich, iaphi0, iaphie
    integer(kind=8) :: iaprng, iaprnl, iarefe, iasupm, iavall, iavalp
    integer(kind=8) :: ibid, iblph, icmp, icog, icol
    integer(kind=8) :: ieqg, ieql, iiblph, ili, inoe, inog
    integer(kind=8) :: inol, iret, isma, j, lgblph, nblph
    integer(kind=8) :: nbnoet, ncmpmx, nddle, nddli, nddlt, nec, nlblph
    integer(kind=8) :: nueqg, nueql
    integer(kind=8), pointer :: desm(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    character(len=24), pointer :: refe(:) => null()
    real(kind=8), pointer :: para_r(:) => null()
    integer(kind=8), pointer :: nueg(:) => null()
    integer(kind=8), pointer :: nuel(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    chno_g = chno_gz
    chno_l = chno_lz
!
    call dismoi('NOM_MAILLA', chno_gz, 'CHAM_NO', repk=mag)
    call jeveuo(chno_g//'.REFE', 'L', vk24=refe)
    nume_equa_g = refe(2) (1:19)
    call nueq_chck(nume_equa_g)
    call dismoi('NOM_GD', chno_gz, 'CHAM_NO', repk=nomgd)
    if (nomgd(1:6) .ne. 'DEPL_R') then
        call utmess('F', 'SOUSTRUC_43', sk=nomgd)
    end if
!
!
!     1- RECUPERATION DU NOM DU MACR_ELEM:
!     ------------------------------------
    call jeveuo(mag//'.NOMACR', 'L', vk8=vnomacr)
    call jenonu(jexnom(mag//'.SUPMAIL', mail), isma)
    if (isma .le. 0) then
        ch8(1) = mail
        ch8(2) = mag
        call utmess('F', 'SOUSTRUC_44', nk=2, valk=ch8)
    end if
    call jeveuo(jexnom(mag//'.SUPMAIL', mail), 'L', iasupm)
    nomacr = vnomacr(isma)
    nul = nomacr
    nume_equa_l = nul//'.NUME'
    call nueq_chck(nume_equa_l)
!
    call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
    call jeveuo(nomacr//'.CONX', 'L', vi=conx)
    call jeveuo(nomacr//'.DESM', 'L', vi=desm)
    nbnoet = desm(2)+desm(8)+desm(9)
    nddle = desm(4)
    nddli = desm(5)
    nddlt = nddle+nddli
!                 '&&SSDEIN.VALP' EST UN VECTEUR DE TRAVAIL :
    call wkvect('&&SSDEIN.VALP', 'V V R', nddlt, iavalp)
!
    call jeveuo(chno_g//'.VALE', 'L', vr=vale)
    call jeveuo(jexnum(nume_equa_g//'.PRNO', 1), 'L', iaprng)
    call jeveuo(nume_equa_g//'.NUEQ', 'L', vi=nueg)
    call jeveuo(nume_equa_l//'.NUEQ', 'L', vi=nuel)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=ncmpmx)
!
!
!     2- ALLOCATION DU CHAM_NO RESULTAT : CHNO_L
!     --------------------------------------
!     .REFE:
    call wkvect(chno_l//'.REFE', 'G V K24', 4, iarefe)
    call jeecra(chno_l//'.REFE', 'DOCU', ibid, 'CHNO')
    zk24(iarefe-1+2) = nul//'.NUME'
!     .VALE:
    call wkvect(chno_l//'.VALE', 'G V R', nddlt, iavall)
!
!
!     4- CALCUL DES VALEURS DE CHNO_L.VALE:
!     ---------------------------------
!
!     4-1- ON RECOPIE CHNO_G.VALE DANS Q_E:
!     ---------------------------------
    do inoe = 1, nbnoet
        inog = zi(iasupm-1+inoe)
        inol = conx(3*(inoe-1)+2)
        ili = conx(3*(inoe-1)+1)
!
        call jeveuo(jexnum(nume_equa_l//'.PRNO', ili), 'L', iaprnl)
!
        nueql = zi(iaprnl-1+(inol-1)*(nec+2)+1)
        iadgl = iaprnl-1+(inol-1)*(nec+2)+3
        ieql = nuel(nueql)
        if (ieql .le. nddli) then
            call utmess('F', 'SOUSTRUC_45')
        end if
!
        nueqg = zi(iaprng-1+(inog-1)*(nec+2)+1)
        iadgg = iaprng-1+(inog-1)*(nec+2)+3
!
        icol = 0
        icog = 0
        do 2, icmp = 1, ncmpmx
            exil = exisdg(zi(iadgl), icmp)
            exig = exisdg(zi(iadgg), icmp)
            if (exil) icol = icol+1
            if (exig) icog = icog+1
            if (exig .and. exil) then
                ieql = nuel(nueql-1+icol)
                ieqg = nueg(nueqg-1+icog)
                zr(iavall-1+ieql) = vale(ieqg)
            end if
2           continue
        end do
!
!
!     4-2- ON CHANGE LE REPERE (ROTATION G->L ) : Q_E  --> Q_E :
!     ----------------------------------------------------------
!
        call ssrone(mag, isma, rota)
!
        if (rota(1:3) .eq. 'OUI') then
            call jeveuo(mag//'.PARA_R', 'L', vr=para_r)
            angl(1) = para_r(14*(isma-1)+4)
            angl(2) = para_r(14*(isma-1)+5)
            angl(3) = para_r(14*(isma-1)+6)
            call matrot(angl, pgl)
            do i = 1, 3
                do j = 1, 3
                    lambda(i, j) = pgl(i, j)
                    lambda(i, j+3) = 0.d0
                    lambda(i+3, j) = 0.d0
                    lambda(i+3, j+3) = pgl(i, j)
                end do
            end do
            call ssvaro(lambda, 'GL', .false._1, 'EXTE', nomacr, &
                        iavall, iavalp)
            do i = 1, nddle
                zr(iavall-1+nddli+i) = zr(iavalp-1+nddli+i)
            end do
            call jedetr('&&SSVARO.IINO')
        end if
!
!
!     4-3  Q_I= (K_II**-1)*F_I :
!     -------------------------
        if (nocas(1:1) .ne. ' ') then
            call jeexin(jexnom(nomacr//'.LICA', nocas), iret)
            if (iret .eq. 0) then
                valk(1) = nocas
                valk(2) = nomacr
                call utmess('A', 'SOUSTRUC_46', nk=2, valk=valk)
            else
                call jeveuo(jexnom(nomacr//'.LICA', nocas), 'L', ialica)
                call jeveuo(jexnom(nomacr//'.LICH', nocas), 'L', ialich)
!
                if (zk8(ialich-1+1) (1:3) .eq. 'NON') then
!
!           -- LE CHARGEMENT N'EST PAS "SUIVEUR" :
                    if (rota(1:3) .eq. 'OUI') then
                        call ssvaro(lambda, 'GL', .false._1, 'TOUS', nomacr, &
                                    ialica, iavalp)
                        call ssvau1(nomacr, iavalp, iavalp)
                        do i = 1, nddli
                            zr(iavall-1+i) = zr(iavalp-1+i)
                        end do
                    else
                        do i = 1, nddli
                            zr(iavall-1+i) = zr(ialica-1+nddlt+i)
                        end do
                    end if
!
                else if (zk8(ialich-1+1) (1:3) .eq. 'OUI') then
!
!           -- LE CHARGEMENT EST "SUIVEUR" :
                    do i = 1, nddli
                        zr(iavall-1+i) = zr(ialica-1+nddlt+i)
                    end do
                else
                    call utmess('F', 'SOUSTRUC_47')
                end if
            end if
        end if
!
!
!
!     4-4  Q_I= Q_I + PHI_IE * Q_E :
!     ------------------------------
        call jelira(nomacr//'.PHI_IE', 'LONMAX', lgblph)
        call jelira(nomacr//'.PHI_IE', 'NMAXOC', nblph)
        nlblph = lgblph/nddli
!
        j = 0
        do iblph = 1, nblph
            call jeveuo(jexnum(nomacr//'.PHI_IE', iblph), 'L', iaphi0)
            do iiblph = 1, nlblph
                j = j+1
                if (j .gt. nddle) goto 13
                iaphie = iaphi0+(iiblph-1)*nddli
                do i = 1, nddli
                    zr(iavall-1+i) = zr(iavall-1+i)-zr(iaphie-1+i)*zr( &
                                     iavall-1+nddli+j)
                end do
            end do
13          continue
            call jelibe(jexnum(nomacr//'.PHI_IE', iblph))
        end do
!
        call jedetr('&&SSDEIN.VALP')
!
        call jedema()
        end subroutine
