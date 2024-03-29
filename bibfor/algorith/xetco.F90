! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine xetco(champ, nomcha, nomch0)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescns.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
!
!
    character(len=19) :: cescoi, cesco0
    character(len=24) :: ligrel
    integer :: iad1, iad2, ibid, ima, ispt, icmp
    integer :: jcelk, jcesv0, jcesvi, jcl0, jcli
    integer :: jco0, jcoi, nbcmp0, nbcmpi, nbma, nbma0, nbsp0
    integer :: nbspi
    integer :: jcesk, jconx1, jconx2, jcnli, jcnsvi, jfiss, neq
    integer :: nuno, jlis, nuno2, ieq
    character(len=8) :: noma, nomo, fiss, kbid
    character(len=19) :: cnscoi, nliseq
    character(len=24) :: champ, nomcha, nomch0
! --------------------------------------------------
!
! --- TRAITEMENT PARTICULIER SI RECOPIE ÉTAT COHÉSIF
! --- ON VA COMPARER DEUX CARTES DE TAILLES:
! ---  1- CELLE DU .XCOH DU RÉSULTAT DONNÉ EN ENTRÉE
! ---  2- CELLE DU .XCO0 DU RÉSULTAT SORTANT
!
! IN  CHAMP  : CHAMP NON NUL DONNE A SNL/ETAT_INIT
!              POUR INITIALISER LES VARIABLES INTERNES
!              IL N'A PAS LA BONNE STRUCTURE CAR
!              IL CORRESPOND A L ANCIENNE FISSURE
! IN  NOMCH0 : CHAMP NUL MS AVEC LA BONNE STRUCTURE
! OUT NOMCHA : CHAMP INITIAL PRODUIT
!              NON NUL AVEC LA BONNE STRUCTURE
!
! --------------------------------------------------
    call jeveuo(nomch0(1:19)//'.CELK', 'L', jcelk)
    ligrel = zk24(jcelk)
    cescoi = '&&NMETL1.CESCOI'
    call celces(champ, 'V', cescoi)
    call jeveuo(cescoi//'.CESD', 'L', jcoi)
    call jeveuo(cescoi//'.CESL', 'L', jcli)
    call jeveuo(cescoi//'.CESV', 'L', jcesvi)
    cesco0 = '&&NMETL1.CESCO0'
    call celces(nomch0, 'V', cesco0)
    call jeveuo(cesco0//'.CESD', 'L', jco0)
    call jeveuo(cesco0//'.CESL', 'L', jcl0)
    call jeveuo(cesco0//'.CESV', 'E', jcesv0)
    call jeveuo(cesco0//'.CESK', 'L', jcesk)
!
!   RECUP MAILLAGE ET INFOS
    noma = zk8(jcesk)
    call jeveuo(noma//'.CONNEX', 'L', jconx1)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
!   TRANSFO EN CHAM_NO_S
    cnscoi = '&&NMETL1.CNSCO0'
    call cescns(cescoi, ' ', 'V', cnscoi, ' ', &
                ibid)
    call jeveuo(cnscoi//'.CNSL', 'L', jcnli)
    call jeveuo(cnscoi//'.CNSV', 'L', jcnsvi)
!
!   RÉCUP DE LA NOUVELLE FISSURE
    nomo = ligrel(1:8)
    call jeveuo(nomo//'.FISS', 'L', jfiss)
    fiss = zk8(jfiss)
!
!   RECUP DE LA LISTE DE RELATIONS INF-SUP
    nliseq = fiss//'.LISEQ'
    call jelira(nliseq, 'LONMAX', neq, kbid)
    neq = neq/2
    call jeveuo(nliseq, 'L', jlis)
!
!   RÉCUP DU NOMBRE DE MAILLES
    nbma = zi(jcoi)
    nbma0 = zi(jco0)
    ASSERT(nbma .eq. nbma0)
!
!   BOUCLE SUR LES MAILLES
    do ima = 1, nbma
!
!       NOMBRE DE NOEUDS SUR LA MAILLE
        nbspi = zi(jcoi+5-1+4*(ima-1)+1)
        nbsp0 = zi(jco0+5-1+4*(ima-1)+1)
!
!       NOMBRE DE COMPOSANTES
        nbcmpi = zi(jcoi+5-1+4*(ima-1)+3)
        nbcmp0 = zi(jco0+5-1+4*(ima-1)+3)
!
!       ADRESSE D ECRITURE
        call cesexi('C', jco0, jcl0, ima, 1, &
                    1, 1, iad1)
        call cesexi('C', jcoi, jcli, ima, 1, &
                    1, 1, iad2)
!
!       BOUCLE SUR LES NOEUDS ET COMPOSANTES
        do ispt = 1, nbsp0
            do icmp = 1, nbcmp0
                call cesexi('C', jcoi, jcli, ima, ispt, &
                            1, icmp, iad1)
                call cesexi('C', jco0, jcl0, ima, ispt, &
                            1, icmp, iad2)
!
!               CALCUL COHESIF : PHASE DE PREDICTION
                if (icmp .eq. 3 .and. iad2 .gt. 0) then
                    zr(jcesv0-1+iad2) = 1.d0
                    goto 214
                end if
!
!               SI LES DEUX SONT ATTRIBUES, ON RECOPIE
                if (iad1 .gt. 0 .and. iad2 .gt. 0) then
                    zr(jcesv0-1+iad2) = zr(jcesvi-1+iad1)
                end if
!
!               LE NOUVEAU EST ATTRIBUE MS PAS L'ANCIEN
!               ON REGARDE SI LE NUM GLOBAL EST ATTRIBUE
                if (iad1 .le. 0 .and. iad2 .gt. 0) then
                    nuno = zi(jconx1-1+zi(jconx2+ima-1)+ispt-1)
!
!                   SI OUI, ON RECOPIE
                    if (zl(jcnli-1+nbcmpi*(nuno-1)+icmp)) then
                        zr(jcesv0-1+iad2) = zr(jcnsvi-1+nbcmpi*(nuno-1)+icmp)
!
!                   SINON, ON MET LA VALEUR "MATERIAU SAIN"
!                   SAUF SI ARETE VITALE LIANT A UN NOEUD ATTRIBUE
                    else
                        if (icmp .eq. 1) zr(jcesv0-1+iad2) = 0.d0
                        if (icmp .eq. 2) zr(jcesv0-1+iad2) = -1.d0
                        do ieq = 1, neq
                            if (zi(jlis-1+2*(ieq-1)+1) .eq. nuno) then
                                nuno2 = zi(jlis-1+2*(ieq-1)+2)
                                if (zl(jcnli-1+nbcmpi*(nuno2-1)+icmp)) zr(jcesv0-1+iad2) = &
                                    zr( &
                                    jcnsvi-1+nbcmpi*(nuno2-1 &
                                                     )+icmp &
                                    )
                            else if (zi(jlis-1+2*(ieq-1)+2) .eq. nuno) then
                                nuno2 = zi(jlis-1+2*(ieq-1)+1)
                                if (zl(jcnli-1+nbcmpi*(nuno2-1)+icmp)) zr(jcesv0-1+iad2) = &
                                    zr( &
                                    jcnsvi-1+nbcmpi*(nuno2-1 &
                                                     )+icmp &
                                    )
                            end if
                        end do
                    end if
                end if
!
!               ON REPART D UN ETAT ELASTIQUE
                if (icmp .eq. 2 .and. iad2 .gt. 0) then
                    zr(jcesv0-1+iad2) = -abs(zr(jcesv0-1+iad2))
                end if
214             continue
            end do
        end do
    end do
!
!   ON RETRANSFORME EN CHAM_ELEM NORMAL
!   LES OPTIONS DE PROLONGEMENT DOIVENT ETRE EN ACCORD AVEC
!   LA ROUTINE XMELE3
    call cescel(cesco0, ligrel, 'XCVBCA_MORTAR', 'PCOHES', 'NON', &
                ibid, 'V', nomcha, 'F', ibid)
    call detrsd('CHAMP_GD', cescoi)
    call detrsd('CHAMP_GD', cesco0)
end subroutine
