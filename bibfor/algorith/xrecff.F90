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
subroutine xrecff(fiss, typfis, chfond, basfon, fonoeu, &
                  lnoff, conf)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: lnoff
    character(len=8) :: fiss, typfis, conf
    character(len=24) :: chfond, basfon, fonoeu
!
!
! person_in_charge: samuel.geniaut at edf.fr
!
!       RECUPERATION DE LA LISTE DES POINTS DU FOND DE FISSURE
!       SUR LEQUEL ON VA EFFECTUER LE CALCUL + BASE LOCALE EN FOND
!       DE FISSURE
!
!  IN  : FISS   : SD_FISS_XFEM OU SD_FOND_FISSURE
!  IN  : TYPFIS : TYPE D'OBJET POUR DECRIRE LE FOND DE FISSURE
!                 'FONDFISS' OU 'FISSURE' OU 'THETA'
!  OUT : CHFOND : FOND DE FISSURE SUR LEQUEL ON FERA LE POST TRAITEMENT
!  OUT : BASFON : BASE LOCALE RELATIVE A CHFOND
!  OUT : FONOEU : NOEUD du FOND DE FISSURE (FEM)
!  OUT : LNOFF  : NOMBRE DE POINTS DU FOND CHFOND
!  OUT : CONF   : CONFIGURATION DE LA FISSURE (FEM)
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: numfon, ibid, idepfi, iarrfi, ifon, ibas, inoeu
!
    integer(kind=8) :: i, j, k, nfonu, jfonu, jbasu, jnoeu
    real(kind=8) :: smax, s, s1, s2, xyz1, xyz2
    character(len=24) :: fontmp, bastmp, noeutmp, valk
    integer(kind=8), pointer :: fondmult(:) => null()
    real(kind=8), pointer :: fondfiss(:) => null()
    real(kind=8), pointer :: basefond(:) => null()
    character(len=8), pointer :: fonoeud(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
!     LISTE DES NOEUDS DU FOND DE FISSURE EN FEM
    if (typfis .eq. 'FONDFISS') then
        call jeveuo(fiss//'.FOND.NOEU', 'L', vk8=fonoeud)
    end if
!
!     LISTE DES POINTS DES FONDS DE FISSURES
    call jeveuo(fiss//'.FONDFISS', 'L', vr=fondfiss)
!
!     LISTE DES FONDS MULTIPLES EN XFEM
    if (typfis .eq. 'FISSURE') then
        call jeveuo(fiss//'.FONDMULT', 'L', vi=fondmult)
    end if
!
!     BASE LOCALE EN FOND DE FISSURE
    if (typfis .eq. 'FONDFISS') then
!       CET OBJET N'EXISTE QUE SI CONFIG_INIT='COLLEE'
        call dismoi('CONFIG_INIT', fiss, 'FOND_FISS', repk=conf)
        if (conf .eq. 'COLLEE') then
            call jeveuo(fiss//'.BASEFOND', 'L', vr=basefond)
        end if
    else
        call jeveuo(fiss//'.BASEFOND', 'L', vr=basefond)
    end if
!
!     ------------------------------------------------------------------
!     TRAITEMENT DU MOT-CLE NUME_FOND :
!     RESTRICTION DU FOND ET DE LA BASE AU NUMERO DU FOND DEMANDE
!     ------------------------------------------------------------------
!
!     NUMERO DU FOND A TRAITER
    if (typfis .eq. 'FISSURE') then
        call getvis('THETA', 'NUME_FOND', iocc=1, scal=numfon, nbret=ibid)
    else
        numfon = 1
    end if
!
!     DETERMINATION DU NOMBRE DE NOEUDS EN FOND DE FISSURE
    if (typfis .eq. 'FISSURE') then
        idepfi = fondmult(2*(numfon-1)+1)
        iarrfi = fondmult(2*(numfon-1)+2)
        lnoff = iarrfi-idepfi+1
    else
        idepfi = 1
        call jelira(fiss//'.FOND.NOEU', 'LONMAX', lnoff)
    end if
!
!     CREATION DE NOEUDS TEMPORAIRES
    noeutmp = '&&XREFF.FONNOEU_TEMP'
    if (typfis .eq. 'FONDFISS') then
        call wkvect(noeutmp, 'V V K8', lnoff, inoeu)
        do i = 1, lnoff
            zk8(inoeu-1+i) = fonoeud(i)
        end do
    end if
!
!     CREATION D'UN FOND TEMPORAIRE RESTREINT AU NUMFON
    fontmp = '&&XREFF.FONFIS_TEMP'
    call wkvect(fontmp, 'V V R', lnoff*4, ifon)
    do i = 1, lnoff
        do j = 1, 4
            zr(ifon-1+4*(i-1)+j) = fondfiss(4*(i+(idepfi-1)-1)+j)
        end do
    end do
!
!     CREATION D'UNE BASE TEMPORAIRE RESTREINTE AU NUMFON
    bastmp = '&&XREFF.BASFON_TEMP'
    if ((typfis .eq. 'FISSURE') .or. (conf .eq. 'COLLEE')) then
        call wkvect(bastmp, 'V V R', lnoff*6, ibas)
        do i = 1, lnoff
            do j = 1, 6
                zr(ibas-1+6*(i-1)+j) = basefond(6*(i+(idepfi-1)-1)+j)
            end do
        end do
    end if
!
!     ------------------------------------------------------------------
!     TRAITEMENT DU MOT-CLE NB_POINT_FOND :
!     CREATION DU NOUVEAU FOND ET DE LA NOUVELLE BASE
!     ------------------------------------------------------------------
!
!     DOIT-ON PRENDRE UNE REPARTITION UNIFORME ?
    call getvis('THETA', 'NB_POINT_FOND', iocc=1, scal=nfonu, nbret=ibid)
    if (ibid .eq. 0) nfonu = 0
!
    if (nfonu .gt. 0) then
!
!       SI OUI : MODIFICATION DE LA LISTE DES POINTS DU FOND
!                ET DE LA BASE
!
        ASSERT(nfonu .ge. 2)
!
!       CREATION DU FOND MODIFIE
        call wkvect(chfond, 'V V R', 4*nfonu, jfonu)
!
!       CREATION DE LA BASE MODIFIEE
        if ((typfis .eq. 'FISSURE') .or. (conf .eq. 'COLLEE')) then
            call wkvect(basfon, 'V V R', 6*nfonu, jbasu)
        end if
!
!       CREATION DES NOEUDS MODIFIES
        if (typfis .eq. 'FONDFISS') then
            call wkvect(fonoeu, 'V V K8', nfonu, jnoeu)
        end if
!
!       1ER ET DERNIER POINTS
        do j = 1, 4
            zr(jfonu-1+4*(1-1)+j) = zr(ifon-1+4*(1-1)+j)
            zr(jfonu-1+4*(nfonu-1)+j) = zr(ifon-1+4*(lnoff-1)+j)
        end do
!
        if ((typfis .eq. 'FISSURE') .or. (conf .eq. 'COLLEE')) then
            do j = 1, 6
                zr(jbasu-1+6*(1-1)+j) = zr(ibas-1+6*(1-1)+j)
                zr(jbasu-1+6*(nfonu-1)+j) = zr(ibas-1+6*(lnoff-1)+j)
            end do
        end if
!
        if (typfis .eq. 'FONDFISS') then
            zk8(jnoeu) = 'XXXX'
            zk8(jnoeu+(nfonu-1)) = 'XXXX'
        end if
!
!       NOUVEAUX POINTS
        smax = zr(ifon-1+4*(lnoff-1)+4)
        do i = 2, nfonu-1
            s = (i-1)*smax/(nfonu-1)
            do k = 1, lnoff
                if (zr(ifon-1+4*(k-1)+4) .gt. s) goto 110
            end do
110         continue
!         ON INTERPOLE LES COORD ENTRE CELLES DU SEGMENT [K-1,K]
            s1 = zr(ifon-1+4*(k-1-1)+4)
            s2 = zr(ifon-1+4*(k-1)+4)
            do j = 1, 3
                xyz1 = zr(ifon-1+4*(k-1-1)+j)
                xyz2 = zr(ifon-1+4*(k-1)+j)
                zr(jfonu-1+4*(i-1)+j) = xyz1+(xyz2-xyz1)*(s-s1)/(s2-s1)
            end do
            if ((typfis .eq. 'FISSURE') .or. (conf .eq. 'COLLEE')) then
                do j = 1, 6
                    xyz1 = zr(ibas-1+6*(k-1-1)+j)
                    xyz2 = zr(ibas-1+6*(k-1)+j)
                    zr(jbasu-1+6*(i-1)+j) = xyz1+(xyz2-xyz1)*(s-s1)/(s2-s1)
                end do
            end if
!
            zr(jfonu-1+4*(i-1)+4) = s
        end do
!
!       CREATION DES NOEUDS MODIFIES
        if (typfis .eq. 'FONDFISS') then
            valk = 'XXXX'
!
!           FONOEU : VALEURS MISES A XXXX
            do i = 1, nfonu
                zk8(jnoeu-1+i) = valk
            end do
        end if
!
!       ON ECRASE LNOFF
        lnoff = nfonu
!
    else
!
!       SI NON : ON RECOPIE TELLE QUELLE LA LISTE DES POINTS DU FOND
        call jedupo(fontmp, 'V', chfond, .false._1)
!
!       ET ON RECOPIE TELLE QUELLE LA BASE
        call jedupo(bastmp, 'V', basfon, .false._1)
!
!       ET ON RECOPIE TELS QUELS LES NOEUDS
        if (typfis .eq. 'FONDFISS') then
            call jedupo(noeutmp, 'V', fonoeu, .false._1)
        end if
    end if
!
!     MENAGE
    call jedetr(fontmp)
    call jedetr(bastmp)
    call jedetr(noeutmp)
!
    call jedema()
end subroutine
