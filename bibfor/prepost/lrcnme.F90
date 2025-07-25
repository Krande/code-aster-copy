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
subroutine lrcnme(chanom, nochmd, nomamd, nomaas, nomgd, &
                  typent, nbcmpv, ncmpva, ncmpvm, iinst, &
                  numpt, numord, inst, crit, prec, &
                  nrofic, codret, base)
!_____________________________________________________________________
!
! person_in_charge: nicolas.sellenet at edf.fr
!     LECTURE D'UN CHAMP AUX NOEUDS - FORMAT MED
!     -    -       -         -               --
!-----------------------------------------------------------------------
!     ENTREES:
!        CHANOM : NOM ASTER DU CHAMP A LIRE
!        NOCHMD : NOM MED DU CHAMP DANS LE FICHIER
!        NOMAMD : NOM MED DU MAILLAGE LIE AU CHAMP A LIRE
!                  SI ' ' : ON SUPPOSE QUE C'EST LE PREMIER MAILLAGE
!                           DU FICHIER
!        TYPENT : TYPE D'ENTITE DU CHAMP
!                (MED_NOEUD=3,MED_MAILLE=0,MED_NOEUD_MAILLE=4)
!        NOMAAS : NOM ASTER DU MAILLAGE
!        NOMGD  : NOM DE LA GRANDEUR ASSOCIEE AU CHAMP
!        NBCMPV : NOMBRE DE COMPOSANTES VOULUES
!                 SI NUL, ON LIT LES COMPOSANTES A NOM IDENTIQUE
!        NCMPVA : LISTE DES COMPOSANTES VOULUES POUR ASTER
!        NCMPVM : LISTE DES COMPOSANTES VOULUES DANS MED
!        IINST  : 1 SI LA DEMANDE EST FAITE SUR UN INSTANT, 0 SINON
!        NUMPT  : NUMERO DE PAS DE TEMPS EVENTUEL
!        NUMORD : NUMERO D'ORDRE EVENTUEL DU CHAMP
!        INST   : INSTANT EVENTUEL
!        CRIT   : CRITERE SUR LA RECHERCHE DU BON INSTANT
!        PREC   : PRECISION SUR LA RECHERCHE DU BON INSTANT
!        NROFIC : NUMERO NROFIC LOGIQUE DU FICHIER MED
!     SORTIES:
!        CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_____________________________________________________________________
!
    use as_med_module, only: as_med_open
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_mfdfdi.h"
#include "asterfort/as_mfdnfc.h"
#include "asterfort/as_mfdnfd.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/cnscno.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lrcame.h"
#include "asterfort/ulisog.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: chanom
    character(len=*) :: ncmpva, ncmpvm
    character(len=8) :: nomaas, nomgd
    character(len=8) :: crit
    character(len=*) :: nochmd, nomamd
!
    integer(kind=8) :: nrofic, typent
    integer(kind=8) :: nbcmpv
    integer(kind=8) :: iinst, numpt, numord
    integer(kind=8) :: codret
    character(len=1), optional, intent(in) :: base
!
    real(kind=8) :: inst
    real(kind=8) :: prec
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='LRCNME')
!
    integer(kind=8) :: iaux, iret, jcmpva, nbcmpa, nbcham, i, nbcmp, jcmp, junit
    integer(kind=8) :: ibid, nseqca, tycha, j
    med_idt :: idfimd
    integer(kind=8) :: jnocmp, ncmprf, ubid
    parameter(ubid=1)
    integer(kind=8) :: unbid(ubid)
    integer(kind=8) :: edlect
    parameter(edlect=0)
!
    character(len=1) :: saux01
    character(len=8) :: saux08, parbid
    character(len=19) :: chamn
    character(len=19) :: chamns, ligbid
    character(len=24) :: optbid
    character(len=64) :: nomcha
    character(len=200) :: nofimd
    character(len=255) :: kfic
    aster_logical :: ttt
    character(len=1) :: bas2
!
!====
! 1. ALLOCATION D'UN CHAM_NO_S  (CHAMNS)
!====
!
    call jemarq()
    if (.not. present(base)) then
        bas2 = 'G'
    else
        bas2 = base
    end if
!
! 1.1. ==> REPERAGE DES CARACTERISTIQUES DE CETTE GRANDEUR
!
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), iaux)
    if (iaux .eq. 0) then
        call utmess('F', 'MED_65')
    end if
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jnocmp)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', ncmprf)
!
! 1.2. ==> ALLOCATION DU CHAM_NO_S
!
!               1234567890123456789
    chamns = '&&      .CNS.MED   '
    chamns(3:8) = nompro
!
    call jeexin(ncmpva, iret)
    if (iret .gt. 0) then
        call jeveuo(ncmpva, 'L', jcmpva)
        call jelira(ncmpva, 'LONMAX', nbcmpa)
        if (nbcmpa .le. ncmprf) then
            do i = 1, nbcmpa
                ttt = .false.
                do j = 1, ncmprf
                    if (zk8(jcmpva+i-1) .eq. zk8(jnocmp+j-1)) ttt = .true.
                end do
                if (.not. ttt) then
                    call utmess('F', 'MED_66')
                end if
            end do
        else
            call utmess('F', 'MED_70')
        end if
!
    else
        call ulisog(nrofic, kfic, saux01)
        if (kfic(1:1) .eq. ' ') then
            call codent(nrofic, 'G', saux08)
            nofimd = 'fort.'//saux08
        else
            nofimd = kfic(1:200)
        end if
        call as_med_open(idfimd, nofimd, edlect, iret)
        call as_mfdnfd(idfimd, nbcham, iret)
        do i = 1, nbcham
            call as_mfdnfc(idfimd, i, nbcmp, iret)
            call wkvect('&&LRCNME.NOMCMP_K16', 'V V K16', nbcmp, jcmp)
            call wkvect('&&LRCNME.UNITCMP', 'V V K16', nbcmp, junit)
            call as_mfdfdi(idfimd, i, nomcha, tycha, zk16(jcmp), &
                           zk16(junit), nseqca, iret)
            if (nomcha .eq. nochmd) then
                ncmprf = nbcmp
                call wkvect('&&LRCNME.NOMCMP_K8', 'V V K8', nbcmp, jnocmp)
                do j = 1, nbcmp
                    zk8(jnocmp+j-1) = zk16(jcmp+j-1) (1:8)
                end do
                call jedetr('&&LRCNME.NOMCMP_K16')
                call jedetr('&&LRCNME.UNITCMP')
                goto 780
            end if
            call jedetr('&&LRCNME.NOMCMP_K16')
            call jedetr('&&LRCNME.UNITCMP')
780         continue
        end do
        call as_mficlo(idfimd, iret)
!
    end if
!
!====
! 2. LECTURE
!====
!
    ligbid = ' '
    optbid = ' '
    parbid = ' '
    call lrcame(nrofic, nochmd, nomamd, nomaas, ligbid, &
                optbid, parbid, 'NOEU', typent, unbid, &
                unbid, unbid, nbcmpv, ncmpva, ncmpvm, &
                iinst, numpt, numord, inst, crit, &
                prec, nomgd, ncmprf, jnocmp, chamns, &
                codret)
!
!====
! 3. TRANSFORMATION DU CHAM_NO_S EN CHAM_NO :
!====
!
    chamn = chanom
!
    call cnscno(chamns, ' ', 'NON', bas2, chamn, &
                'F', ibid)
!
    call detrsd('CHAM_NO_S', chamns)
!
!====
! 4. BILAN
!====
!
    if (codret .ne. 0) then
        call utmess('A', 'MED_55', sk=chamn)
    end if
    call jedetr('&&LRCNME.NOMCMP_K8')
    call jedema()
!
end subroutine
