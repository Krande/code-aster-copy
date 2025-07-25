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
subroutine xpocmp(elrefp, cns1, ima, n, jconx1, &
                  jconx2, ndim, nfh, nfe, ddlc, &
                  nbcmp, cmp, lmeca, pre1)
!
! person_in_charge: samuel.geniaut at edf.fr
!
! aslint: disable=W1306
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: ndim, nfh, nfe, ima, n, jconx1, jconx2, nbcmp, cmp(*)
    integer(kind=8) :: ddlc
    aster_logical :: lmeca, pre1, press, press1, pref, lagf
    character(len=8) :: elrefp
    character(len=19) :: cns1
!
!   DETERMINER LES COMPOSANTES ACTIVES DU CHAMP DE DEPLACEMENT
!
!   IN
!     CNS1   : CHAMP_NO_S DU DEPLACEMENT EN ENTREE
!     IMA    : NUMERO DE MAILLE COURANTE PARENT
!     N      : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
!     JCONX1 : ADRESSE DE LA CONNECTIVITE DU MAILLAGE SAIN
!     JCONX2 : LONGUEUR CUMULEE DE LA CONNECTIVITE DU MAILLAGE SAIN
!     NDIM   : DIMENSION DU MAILLAGE
!     NBCMP  : NOMBRE DE COMPOSANTES DU CHAMP_NO DE DEPL1
!     LMECA  : VRAI DANS LE CAS MECANIQUE (SINON CAS THERMIQUE)
!
!   OUT
!     NFH    : NOMBRE DE FONCTIONS HEAVISIDE (PAR NOEUD)
!     NFE    : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT (1 A 4)
!     CMP    : POSITION DES DDLS DE DEPL X-FEM DANS LE CHAMP_NO DE DEPL1
!
!
    integer(kind=8) :: jcnsl1, i, j, k, ino, icmp, ndc, ipos, nnos, ibid
    aster_logical :: exist(n, nbcmp), contas
    character(len=8) :: nomcmp, k8bid
    character(len=8), pointer :: cnsc(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     COMPOSANTES DU CHAMP DE DEPLACEMENT 1
    call jeveuo(cns1//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(cns1//'.CNSL', 'L', jcnsl1)
!
    do j = 1, n
!       INO : NUMÉRO DU NOEUD DANS MALINI
        ino = zi(jconx1-1+zi(jconx2+ima-1)+j-1)
        do icmp = 1, nbcmp
            exist(j, icmp) = zl(jcnsl1-1+(ino-1)*nbcmp+icmp)
        end do
    end do
!
!     ON REGARDE LES COMPOSANTES ACTIVES EN CHAQUE NOEUD
    ipos = 0
    ndc = 0
    nfh = 0
    nfe = 0
    ddlc = 0
    contas = .true.
    press = .true.
    press1 = .true.
    pref = .true.
    lagf = .true.
    call elelin(1, elrefp, k8bid, ibid, nnos)
!
    do i = 1, nbcmp
        nomcmp = cnsc(i)
!
        if (nomcmp(1:4) .eq. 'LAGS' .or. nomcmp(1:4) .eq. 'LAG2' .or. nomcmp(1:4) .eq. &
            'LAG3' .or. nomcmp(1:4) .eq. 'LAG4' .or. nomcmp(1:2) .eq. 'D1' .or. nomcmp(1:2) &
            .eq. 'V1' .or. nomcmp(1:2) .eq. 'D2' .or. nomcmp(1:2) .eq. 'V2' .or. &
            nomcmp(1:2) .eq. 'D3' .or. nomcmp(1:2) .eq. 'V3') then
            do k = 1, nnos
                if (.not. exist(k, i)) contas = .false.
            end do
            if (contas) goto 1
        end if
!
        if (nomcmp(1:4) .eq. 'PRE1') then
            do k = 1, nnos
                if (.not. exist(k, i)) press = .false.
            end do
            if (press) goto 1
        end if
!
        if (nomcmp(1:6) .eq. 'H1PRE1' .or. nomcmp(1:6) .eq. 'H2PRE1' .or. nomcmp(1:6) .eq. &
            'H3PRE1') then
            do k = 1, nnos
                if (.not. exist(k, i)) press1 = .false.
            end do
            if (press1) goto 1
        end if
!
        if (nomcmp(1:7) .eq. 'PRE_FLU' .or. nomcmp(1:7) .eq. 'PR2_FLU' .or. nomcmp(1:7) &
            .eq. 'PR3_FLU') then
            do k = 1, nnos
                if (.not. exist(k, i)) pref = .false.
            end do
            if (pref) goto 1
        end if
!
        if (nomcmp(1:6) .eq. 'LAG_FL' .or. nomcmp(1:6) .eq. 'LA2_FL' .or. nomcmp(1:6) .eq. &
            'LA3_FL') then
            do k = 1, nnos
                if (.not. exist(k, i)) lagf = .false.
            end do
            if (lagf) goto 1
        end if
!
        do j = 1, n
            if (.not. exist(j, i)) goto 21
        end do
!
1       continue
!
        if (nomcmp(1:2) .eq. 'DX' .or. nomcmp(1:2) .eq. 'DY' .or. nomcmp(1:2) .eq. 'DZ' &
            .or. nomcmp(1:1) .eq. 'T') then
            ipos = ipos+1
            ndc = ndc+1
            cmp(ipos) = i
        end if
        if (pre1) then
            if (nomcmp(1:4) .eq. 'PRE1') then
                ipos = ipos+1
                cmp(ipos) = i
            end if
        end if
        if (nomcmp(1:1) .eq. 'H' .and. nomcmp(3:3) .ne. 'P') then
            ipos = ipos+1
            nfh = nfh+1
            cmp(ipos) = i
        end if
        if (pre1) then
            if (nomcmp(1:6) .eq. 'H1PRE1' .or. nomcmp(1:6) .eq. 'H2PRE1' .or. nomcmp(1:6) &
                .eq. 'H3PRE1') then
                ipos = ipos+1
                cmp(ipos) = i
            end if
        end if
        if (nomcmp(1:2) .eq. 'K1' .or. nomcmp(1:2) .eq. 'K2' .or. nomcmp(1:2) .eq. 'K3') then
            ipos = ipos+1
            nfe = nfe+1
            cmp(ipos) = i
        end if
        if (nomcmp(1:2) .eq. 'E1') then
            ASSERT(.not. lmeca)
            nfe = nfe+1
            ipos = ipos+1
            cmp(ipos) = i
        end if
        if (.not. pre1) then
            if (nomcmp(1:3) .eq. 'LAG') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
        end if
        if (pre1) then
            if (nomcmp(1:7) .eq. 'PRE_FLU' .or. nomcmp(1:7) .eq. 'PR2_FLU' .or. nomcmp(1:7) &
                .eq. 'PR3_FLU') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
            if (nomcmp(1:6) .eq. 'LAG_FL' .or. nomcmp(1:6) .eq. 'LA2_FL' .or. nomcmp(1:6) &
                .eq. 'LA3_FL') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
            if (nomcmp(1:4) .eq. 'LAGS' .or. nomcmp(1:4) .eq. 'LAG2' .or. nomcmp(1:4) .eq. &
                'LAG3') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
            if (nomcmp(1:2) .eq. 'D1' .or. nomcmp(1:2) .eq. 'D2' .or. nomcmp(1:2) .eq. 'D3') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
            if (nomcmp(1:2) .eq. 'V1' .or. nomcmp(1:2) .eq. 'V2' .or. nomcmp(1:2) .eq. 'V3') then
                ipos = ipos+1
                ddlc = ddlc+1
                cmp(ipos) = i
            end if
        end if
!
21      continue
    end do
!
    if (lmeca) then
!       CAS DE LA MECANIQUE
        nfe = nfe/ndim
        nfh = nfh/ndim
        ASSERT(ndim .eq. ndc)
    else
!       CAS DE LA THERMIQUE
        ASSERT(ndc .eq. 1)
    end if
!
    call jedema()
end subroutine
