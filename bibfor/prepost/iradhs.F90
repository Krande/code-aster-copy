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
subroutine iradhs(versio)
!
    implicit none
!
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/inistb.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utidea.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8), intent(in) :: versio
!
! --------------------------------------------------------------------------------------------------
!
!                BUT: TRAITER LES "ADHERENCES SUPERTAB":
!     OBJETS JEVEUX CREES:
!        &&IRADHS.CODEGRA : TABLEAU D'ENTIERS DIMENSIONNE
!                 AU NOMBRE DE TYPE_MAILLES DONNANT LE CODE GRAPHIQUE
!                 DE CHAQUE TYPE_MAILLE POUR SUPERTAB.
!        &&IRADHS.CODEPHY : TABLEAU D'ENTIERS DIMENSIONNE
!                 AU NOMBRE DE TYPE_ELEM DONNANT LE CODE PHYSIQUE
!                 DE CHAQUE TYPE_ELEMENT POUR SUPERTAB.
!        &&IRADHS.CODEPHD : TABLEAU D'ENTIERS DIMENSIONNE
!                 AU NOMBRE DE TYPE_MAILLE DONNANT LE CODE PHYSIQUE
!                 PAR DEFAUT DE CHAQUE TYPE_MAILLE POUR SUPERTAB.
!        &&IRADHS.PERMUTA : TABLEAU D'ENTIERS DIMENSIONNE
!                 AU NOMBRE DE TYPE_MAILLES*NOMBRE MAXI DE NOEUDS/MAILLE
!                 DONNANT LES PERMUTATIONS EVENTUELLES DANS L'ORDRE DES
!                 NOEUDS DE CHAQUE TYPE_MAILLE, ENTRE ASTER ET SUPERTAB.
!                 LE NOMBRE MAXI DE NOEUDS/MAILLE EST STOCKE AU BOUT DU
!                 VECTEUR
!
!   IN:  VERSIO = VERSION D'IDEAS 4 OU 5(DEFAUT)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=2), parameter :: axdpcp(4) = (/'AX', 'DP', 'CP', 'PL'/)
    character(len=5) :: phe(2), mot
    character(len=8) :: nommai
    character(len=16) :: nomele
    integer(kind=8) ::  iax, iel, ima, imper, ino, inos
    integer(kind=8) :: iphe, iret1, iret2, iret3, iret4, iret5, iret6
    integer(kind=8) :: isu, itel, jcodd, jpefsu
    integer(kind=8) :: jpermu, jpersu, nbn, nbtyel, nbtyma
    integer(kind=8) :: nbtyms
    integer(kind=8), parameter :: maxnod = 32
    integer(kind=8), parameter :: maxfa = 6
    character(len=8) :: nomtm
    integer(kind=8) :: icas
    integer(kind=8), pointer :: codephy(:) => null()
    character(len=8), pointer :: typema(:) => null()
    integer(kind=8), pointer :: codegra(:) => null()
    integer(kind=8), pointer :: nbno(:) => null()
!
    character(len=8) :: nomail(MT_NTYMAX)
    integer(kind=8) :: limail(MT_NTYMAX), indic(MT_NTYMAX), indicf(MT_NTYMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call jeexin('&&IRADHS.PERMSUP', iret1)
    call jeexin('&&IRADHS.PERFSUP', iret2)
    call jeexin('&&IRADHS.CODEGRA', iret3)
    call jeexin('&&IRADHS.PERMUTA', iret4)
    call jeexin('&&IRADHS.CODEPHY', iret5)
    call jeexin('&&IRADHS.CODEPHD', iret6)
    if (iret1*iret2*iret3*iret4*iret5*iret6 .ne. 0) goto 999
!
    if (iret1 .eq. 0) then
        call wkvect('&&IRADHS.PERMSUP', 'V V I', maxnod*MT_NTYMAX, jpersu)
    end if
    call jeveuo('&&IRADHS.PERMSUP', 'E', jpersu)
    if (iret2 .eq. 0) then
        call wkvect('&&IRADHS.PERFSUP', 'V V I', maxfa*MT_NTYMAX, jpefsu)
    end if
    call jeveuo('&&IRADHS.PERFSUP', 'E', jpefsu)
    call inistb(maxnod, nbtyms, nomail, indic, zi(jpersu), &
                limail, indicf, zi(jpefsu), maxfa)
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtyma)
    if (iret3 .eq. 0) then
        call jecreo('&&IRADHS.CODEGRA', 'V V I')
        call jeecra('&&IRADHS.CODEGRA', 'LONMAX', nbtyma)
    end if
    call jeveuo('&&IRADHS.CODEGRA', 'E', vi=codegra)
    do ima = 1, nbtyma
        call jenuno(jexnum('&CATA.TM.NOMTM', ima), nomtm)
        if (nomtm .eq. 'HEXA27') nomtm = 'HEXA20'
        if (nomtm .eq. 'PENTA18') nomtm = 'PENTA15'
        if (nomtm .eq. 'TRIA7') nomtm = 'TRIA6'
        if (nomtm .eq. 'QUAD9') nomtm = 'QUAD8'
        if (nomtm .eq. 'SEG4') nomtm = 'SEG2'
        do isu = 1, nbtyms
            if (nomtm .eq. nomail(isu)) then
                codegra(ima) = isu
                goto 1
            end if
        end do
1       continue
    end do
    if (iret4 .eq. 0) then
        call wkvect('&&IRADHS.PERMUTA', 'V V I', maxnod*MT_NTYMAX+1, jpermu)
        zi(jpermu-1+maxnod*MT_NTYMAX+1) = maxnod
    end if
    if (iret6 .eq. 0) then
        call wkvect('&&IRADHS.CODEPHD', 'V V I', nbtyma, jcodd)
    end if
    call jeveuo('&&IRADHS.PERMUTA', 'E', jpermu)
    call jeveuo('&CATA.TM.NBNO', 'L', vi=nbno)
    do ima = 1, nbtyma
        nbn = nbno(ima)
        isu = codegra(ima)
        if (isu .eq. 0) then
            icas = 0
        else
            icas = indic(isu)
        end if
!
        if (icas .lt. 0) then
            do ino = 1, nbn
                zi(jpermu-1+maxnod*(ima-1)+ino) = 0
            end do
        else if (icas .eq. 0) then
            do ino = 1, nbn
                zi(jpermu-1+maxnod*(ima-1)+ino) = ino
            end do
        else
            do ino = 1, nbn
                do inos = 1, nbn
                    imper = zi(jpersu-1+maxnod*(isu-1)+inos)
                    if (ino .eq. imper) then
                        zi(jpermu-1+maxnod*(ima-1)+ino) = inos
                        goto 43
                    end if
                end do
43              continue
            end do
        end if
        call jenuno(jexnum('&CATA.TM.NOMTM', ima), nommai)
        call utidea(nommai, zi(jcodd-1+ima), versio)
    end do
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nbtyel)
    if (iret5 .eq. 0) then
        call jecreo('&&IRADHS.CODEPHY', 'V V I')
        call jeecra('&&IRADHS.CODEPHY', 'LONMAX', nbtyel)
        call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typema)
    end if
    call jeveuo('&&IRADHS.CODEPHY', 'E', vi=codephy)
    do itel = 1, nbtyel
        call jenuno(jexnum('&CATA.TE.NOMTE', itel), nomele)
        call utidea(typema(itel), codephy(itel), versio)
    end do
    call jenonu(jexnom('&CATA.TE.NOMTE', 'MEDKQU4'), iel)
    if (iel .ne. 0) codephy(iel) = 94
    call jenonu(jexnom('&CATA.TE.NOMTE', 'MEDKTR3'), iel)
    if (iel .ne. 0) codephy(iel) = 91
    call jenonu(jexnom('&CATA.TE.NOMTE', 'MEDSQU4'), iel)
    if (iel .ne. 0) codephy(iel) = 94
    call jenonu(jexnom('&CATA.TE.NOMTE', 'MEDSTR3'), iel)
    if (iel .ne. 0) codephy(iel) = 91
    call jenonu(jexnom('&CATA.TE.NOMTE', 'MEQ4QU4'), iel)
    if (iel .ne. 0) codephy(iel) = 94
    phe(1) = 'MECA_'
    phe(2) = 'THER_'
    do iphe = 1, 2
        do iax = 1, 4
            mot = phe(iphe) (1:2)//axdpcp(iax)
            call jenonu(jexnom('&CATA.TE.NOMTE', mot(1:4)//'QU4'), iel)
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'AX') codephy(iel) = 84
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'CP') codephy(iel) = 44
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'DP') codephy(iel) = 54
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'PL') codephy(iel) = 44
            call jenonu(jexnom('&CATA.TE.NOMTE', mot(1:4)//'QU8'), iel)
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'AX') codephy(iel) = 85
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'CP') codephy(iel) = 45
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'DP') codephy(iel) = 55
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'PL') codephy(iel) = 45
            call jenonu(jexnom('&CATA.TE.NOMTE', mot(1:4)//'TR3'), iel)
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'AX') codephy(iel) = 81
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'CP') codephy(iel) = 41
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'DP') codephy(iel) = 51
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'PL') codephy(iel) = 41
            call jenonu(jexnom('&CATA.TE.NOMTE', mot(1:4)//'TR6'), iel)
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'AX') codephy(iel) = 82
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'CP') codephy(iel) = 42
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'DP') codephy(iel) = 52
            if (iel .ne. 0 .and. axdpcp(iax) .eq. 'PL') codephy(iel) = 42
            call jenonu(jexnom('&CATA.TE.NOMTE', mot(1:4)//'SE2'), iel)
            if (iel .ne. 0) codephy(iel) = 21
        end do
    end do
    call jedetr('&&IRADHS.PERMSUP')
    call jedetr('&&IRADHS.PERFSUP')
999 continue
    call jedema()
end subroutine
