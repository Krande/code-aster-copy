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
subroutine rvgnoe(mcf, iocc, nmaila, nlstnd, nbtrou, &
                  linoeu)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/oreino.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: iocc, nbtrou, linoeu(*)
    character(len=*) :: mcf
    character(len=8) :: nmaila
    character(len=24) :: nlstnd
!     SAISIE DES NOEUDS DE L' OCCURENCE IOCC DE ACTION
!     CAS OU LE LIEU EST UNE LISTE DE GROUP_NO ET/OU NOEUD
!     ------------------------------------------------------------------
! IN   IOCC   : I : NUMERO DE L' OCCURENCE TRAITEE
! IN   NMAILA : K : NOM DU MAILLAGE CONTENANT LES GROUPES ET LES NOEUDS
! JXIN NLSTND : K : NOM OJB S V I <-- NUMERO DES NOEUDS
!     ------------------------------------------------------------------
!     CONSTRUCTION DE LA LISTE NLSTND :
!         LA LISTE ARGUMENT DE NOEUD
!         LES NOEUD DES GROUPES DE NOEUD DANS L' ORDRE DE LA LISTE
!         ARGUMENT DE GROUP_NO
!         SI 2 NOEUD CONSECUTIF DANS CETTE CONSTRUCTION SONT IDENTIQUES
!         ON N' EN GARDE QU' UN
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbrgpn, nbneud, aneud, agrpn, alndtp, alstnd, agneud
    integer(kind=8) :: i, j, k, libre, numnd, nbtnd, n1, nbn, iret, iera
    integer(kind=8) :: asgtu, i1, i2, ny, ier
    real(kind=8) :: vecty(3), tole
    character(len=8) :: courbe, crit
    character(len=24) :: nomgrn
    character(len=15) :: nrepnd
    character(len=17) :: nrepgn
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: list_n(:) => null()
    aster_logical lnomnoe
!
!==================== CORPS DE LA ROUTINE =============================
!
    call jemarq()
!
    nbtnd = 0
    nrepgn = nmaila//'.GROUPENO'
    nrepnd = nmaila//'.NOMNOE'
    call jeexin(nrepnd, ier)
    lnomnoe = .false.
    if (ier .ne. 0) then
        lnomnoe = .true.
    end if
    libre = 1
!   INDICATEUR DE COMMANDE POUR OREINO: 2-POST_RELEVE_T/PRECISION
    iera = 2
!
! --- RECUPERATION DES ENTITES
!
    call getvem(nmaila, 'GROUP_NO', mcf, 'GROUP_NO', iocc, &
                0, zk24(1), nbrgpn)
    call getvem(nmaila, 'NOEUD', mcf, 'NOEUD', iocc, &
                0, zk8, nbneud)
    nbrgpn = -nbrgpn
    nbneud = -nbneud
    if (nbrgpn .ne. 0) then
        call wkvect('&OP0051.NOM.GRP.ND', 'V V K24', nbrgpn, agrpn)
        call getvem(nmaila, 'GROUP_NO', mcf, 'GROUP_NO', iocc, &
                    nbrgpn, zk24(agrpn), n1)
        do i = 1, nbrgpn, 1
            call jeexin(jexnom(nrepgn, zk24(agrpn+i-1)), ier)
            if (ier .ne. 0) then
                call jelira(jexnom(nrepgn, zk24(agrpn+i-1)), 'LONUTI', n1)
                nbtnd = nbtnd+n1
            end if
        end do
    end if
    if (nbneud .ne. 0) then
        call wkvect('&OP0051.NOM.NOEUD', 'V V K8', nbneud, aneud)
        call getvem(nmaila, 'NOEUD', mcf, 'NOEUD', iocc, &
                    nbneud, zk8(aneud), n1)
        nbtnd = nbtnd+nbneud
    end if
!
    if (nbtnd .ne. 0) then
        call wkvect('&OP0051.LIST.ND.TEMP', 'V V I', nbtnd, alndtp)
        do i = 1, nbtnd, 1
            zi(alndtp+i-1) = 0
        end do
!
        if (nbneud .ne. 0) then
            do i = 1, nbneud, 1
                if (lnomnoe) then
                    call jenonu(jexnom(nrepnd, zk8(aneud+i-1)), numnd)
                else
                    numnd = char8_to_int(zk8(aneud+i-1))
                end if
                zi(alndtp+i-1) = numnd
            end do
        end if
        libre = nbneud+1
!
        if (nbrgpn .ne. 0) then
            do i = 1, nbrgpn, 1
                nomgrn = zk24(agrpn+i-1)
                call jeexin(jexnom(nrepgn, nomgrn), ier)
                if (ier .ne. 0) then
                    call jelira(jexnom(nrepgn, nomgrn), 'LONMAX', nbn)
                    call jeveuo(jexnom(nrepgn, nomgrn), 'L', agneud)
                    do j = 1, nbn, 1
                        zi(alndtp+libre-1+j-1) = zi(agneud+j-1)
                    end do
                    libre = libre+nbn
                end if
            end do
        end if
!
        AS_ALLOCATE(vi=list_n, size=libre-1)
!
        nbtnd = 0
        if (nbtrou .eq. 0) then
            list_n(1) = zi(alndtp)
            nbtnd = nbtnd+1
            do i = 1, libre-2, 1
                do j = 1, nbtnd, 1
                    if (zi(alndtp+i) .eq. list_n(j)) goto 250
                end do
                nbtnd = nbtnd+1
                list_n(nbtnd) = zi(alndtp+i)
250             continue
            end do
            libre = nbtnd+1
        else
            do i = 1, nbtrou, 1
                do j = 1, libre-1, 1
                    if (linoeu(i) .eq. zi(alndtp+j-1)) then
                        nbtnd = nbtnd+1
                        goto 252
                    end if
                end do
252             continue
            end do
        end if
!
        if (nbtnd .eq. 0) then
            call utmess('F', 'POSTRELE_64')
        end if
        call wkvect(nlstnd, 'V V I', nbtnd, alstnd)
!
        if (nbtrou .eq. 0) then
            do i = 1, nbtnd, 1
                zi(alstnd+i-1) = list_n(i)
            end do
        else
            nbtnd = libre-1
            libre = 1
            do i = 1, nbtnd, 1
                numnd = zi(alndtp+i-1)
                do j = 1, nbtrou, 1
                    if (linoeu(j) .eq. numnd) then
                        do k = 1, libre-1, 1
                            if (numnd .eq. zi(alstnd+k-1)) goto 302
                        end do
                        zi(alstnd+libre-1) = numnd
                        libre = libre+1
                    end if
                end do
302             continue
            end do
        end if
!
        AS_DEALLOCATE(vi=list_n)
!
!     --- CAS PARTICULIER
!
        call getvr8('ACTION', 'VECT_Y', iocc=iocc, nbval=3, vect=vecty, &
                    nbret=ny)
        if (ny .ne. 0) then
!           VERIFICATIONS PRELIMINAIRES
            if ((nbneud .ge. 2 .and. nbrgpn .eq. 0) .or. (nbneud .eq. 0 .and. nbrgpn .eq. 1)) then
                if (nbrgpn .eq. 1) then
                    nomgrn = zk24(agrpn+1-1)
                    call jeexin(jexnom(nrepgn, nomgrn), ier)
                    if (ier .ne. 0) then
                        call jelira(jexnom(nrepgn, nomgrn), 'LONMAX', nbn)
                        if (nbn .lt. 2) then
                            call utmess('F', 'POSTRELE_21')
                        end if
                    end if
                end if
            else
                call utmess('F', 'POSTRELE_22')
            end if
            call jeexin('&&YAPAS '//'S1   '//'.DESC', n1)
            if (n1 .ne. 0) call jedetr('&&YAPAS '//'S1   '//'.DESC')
            courbe = '&&YAPAS'
            call wkvect(courbe//'S1   '//'.DESC', 'V V R', 6, asgtu)
!           ORIGINE
            i1 = zi(alstnd-1+1)
!           EXTREMITE
            i2 = zi(alstnd-1+libre-1)
            call jeveuo(nmaila//'.COORDO    .VALE', 'L', vr=vale)
!           TOLERANCE
            call getvtx(mcf, 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
            call getvr8(mcf, 'PRECISION', iocc=iocc, scal=tole, nbret=n1)
!           VERIFICATION QUE LES POINTS SONT ALIGNES
            call oreino(nmaila, zi(alstnd), libre-1, i1, i2, &
                        vale, crit, tole, iera, iret)
            if (iret .ne. 0) then
                call utmess('F', 'POSTRELE_60')
            end if
            zr(asgtu-1+1) = vale(3*(i1-1)+1)
            zr(asgtu-1+2) = vale(3*(i1-1)+2)
            zr(asgtu-1+3) = vale(3*(i1-1)+3)
            zr(asgtu-1+4) = vale(3*(i2-1)+1)
            zr(asgtu-1+5) = vale(3*(i2-1)+2)
            zr(asgtu-1+6) = vale(3*(i2-1)+3)
        end if
!
    end if
    call jeexin('&OP0051.NOM.NOEUD', n1)
    if (n1 .ne. 0) call jedetr('&OP0051.NOM.NOEUD')
    call jeexin('&OP0051.NOM.GRP.ND', n1)
    if (n1 .ne. 0) call jedetr('&OP0051.NOM.GRP.ND')
    call jedetr('&OP0051.LIST.ND.TEMP')
!
    call jedema()
end subroutine
