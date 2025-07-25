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

subroutine deelpo(caelem, noma, numail, phie)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8prem.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
    integer(kind=8) :: numail
    real(kind=8) :: phie
    character(len=8) :: caelem, noma
! RECUPERATION DU DIAMETRE EXTERIEUR D'UN ELEMENT DE POUTRE
!-----------------------------------------------------------------------
!  IN : CAELEM : NOM DU CONCEPT DE TYPE CARA_ELEM
!  IN : NOMA   : NOM DU CONCEPT DE TYPE MAILLAGE
!  IN : NUMAIL : NUMERO DE LA MAILLE CORRESPONDANTE
!  OUT: PHIE   : DIAMETRE EXTERIEUR SUR L'ELEMENT
!-----------------------------------------------------------------------
!
!
    integer(kind=8) :: ias, iasbon, iasedi, iasmax, icad, icav, icmp
    integer(kind=8) :: icode, iglma, igrand, ima, inomcp, irang1
    integer(kind=8) :: irang2, iranv1, iranv2, iret, nbcmp, nbma
    integer(kind=8) :: nuenti, numa, nbec
!
!
    real(kind=8) :: difr, r1, r2, tolr
!
    character(len=8) :: nomcmp(2), nomail
    character(len=19) :: carte
    character(len=24) :: cadesc, cavale, calima, gpmama
!
    data nomcmp/'R1      ', 'R2      '/
!
!-----------------------------------------------------------------------
    call jemarq()
!
!-----1.ACCES AUX OBJETS UTILES
!
    gpmama = noma//'.GROUPEMA'
!
    carte = caelem//'.CARGEOPO'
    cadesc = carte//'.DESC'
    cavale = carte//'.VALE'
    calima = carte//'.LIMA'
    call jeexin(cadesc, iret)
    if (iret .eq. 0) then
        call utmess('F', 'ALGELINE_33')
    end if
    call jeveuo(cadesc, 'L', icad)
    call jeveuo(cavale, 'L', icav)
!
!     DETERMINATION DU NUMERO D'ASSOCIATION CORRESPONDANT DANS LA CARTE
    iasmax = zi(icad+1)
    iasedi = zi(icad+2)
    iasbon = 0
!
    do ias = 1, iasedi
        icode = zi(icad+3+2*(ias-1))
        nuenti = zi(icad+3+2*(ias-1)+1)
!        SI AFFECATION SUR UN GROUPE DE MAILLE
        if (icode .eq. 2) then
            call jeveuo(jexnum(gpmama, nuenti), 'L', iglma)
            call jelira(jexnum(gpmama, nuenti), 'LONMAX', nbma)
!        SI AFFECATION SUR UNE LISTE DE MAILLE
        else if (icode .eq. 3) then
            call jeveuo(jexnum(calima, nuenti), 'L', iglma)
            call jelira(jexnum(calima, nuenti), 'LONMAX', nbma)
        end if
!        RECHERCHE DE LA MAILLE
        do ima = 1, nbma
            numa = zi(iglma+ima-1)
            if (numa .eq. numail) then
                iasbon = ias
                goto 40
            end if
        end do
!
    end do
!
    if (iasbon .eq. 0) then
        nomail = int_to_char8(numail)
        call utmess('F', 'ALGELINE_34', sk=nomail)
    end if
!
40  continue
!
!     EXTRACTION DES RAYONS EXTERIEURS AUX DEUX EXTREMITES DE L'ELEMENT
!       SI LE RAYON EXTERIEUR EST CONSTANT SUR L'ELEMENT, ON DEDUIT
!       LE DIAMETRE EXTERIEUR
!
    igrand = zi(icad)
    call jelira(jexnum('&CATA.GD.NOMCMP', igrand), 'LONMAX', nbcmp)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', 'CAGEPO_R'), 'L', inomcp)
!     NOMBRE D'ENTIERS CODES DANS LA CARTE
    call dismoi('NB_EC', 'CAGEPO_R', 'GRANDEUR', repi=nbec)
    irang1 = indik8(zk8(inomcp), nomcmp(1), 1, nbcmp)
    irang2 = indik8(zk8(inomcp), nomcmp(2), 1, nbcmp)
    if (irang1 .eq. 0 .or. irang2 .eq. 0) then
        call utmess('F', 'ALGELINE_35')
    end if
    icode = zi(icad-1+3+2*iasmax+nbec*(iasbon-1)+1)
    iranv1 = 0
    do icmp = 1, irang1
        if (exisdg([icode], icmp)) iranv1 = iranv1+1
    end do
    iranv2 = 0
    do icmp = 1, irang2
        if (exisdg([icode], icmp)) iranv2 = iranv2+1
    end do
    if (iranv1 .eq. 0 .or. iranv2 .eq. 0) then
        call utmess('F', 'ALGELINE_36')
    end if
!
    r1 = zr(icav+nbcmp*(iasbon-1)+iranv1-1)
    r2 = zr(icav+nbcmp*(iasbon-1)+iranv2-1)
    if (r1 .eq. 0.d0 .or. r2 .eq. 0.d0) then
        call utmess('F', 'ALGELINE_37')
    end if
    tolr = r8prem()
    difr = dble(abs(r1-r2))
    if (difr .gt. r1*tolr) then
        call utmess('F', 'ALGELINE_38')
    end if
!
    phie = 2.d0*r1
!
    call jedema()
!
end subroutine
