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
subroutine lrceme(chanom, nochmd, typech, nomamd, nomaas, &
                  nommod, nomgd, typent, nbcmpv, ncmpva, &
                  ncmpvm, prolz, iinst, numpt, numord, &
                  inst, crit, prec, nrofic, option, &
                  param, nbpgma, nbpgmm, nbspmm, codret)
! person_in_charge: nicolas.sellenet at edf.fr
!     LECTURE D'UN CHAMP AUX ELEMENTS - FORMAT MED
!     -    -       -         -               --
!-----------------------------------------------------------------------
!     ENTREES:
!        CHANOM : NOM ASTER DU CHAMP A LIRE
!        NOCHMD : NOM MED DU CHAMP DANS LE FICHIER
!        TYPECH : TYPE DE CHAMP AUX ELEMENTS : ELEM/ELGA/ELNO
!        TYPENT : TYPE D'ENTITE DU CHAMP
!                (MED_NOEUD,MED_MAILLE,MED_NOEUD_MAILLE)
!        NOMAMD : NOM MED DU MAILLAGE LIE AU CHAMP A LIRE
!                  SI ' ' : ON SUPPOSE QUE C'EST LE PREMIER MAILLAGE
!                           DU FICHIER
!        NOMAAS : NOM ASTER DU MAILLAGE
!        NOMMOD : NOM ASTER DU MODELE NECESSAIRE POUR LIGREL
!        NOMGD  : NOM DE LA GRANDEUR ASSOCIEE AU CHAMP
!        NBCMPV : NOMBRE DE COMPOSANTES VOULUES
!                 SI NUL, ON LIT LES COMPOSANTES A NOM IDENTIQUE
!        NCMPVA : LISTE DES COMPOSANTES VOULUES POUR ASTER
!        NCMPVM : LISTE DES COMPOSANTES VOULUES DANS MED
!        PROLZ  : VALEUR DE PROL_ZERO ('OUI' OU 'NAN')
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
! aslint: disable=W1504
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
#include "asterfort/cescar.h"
#include "asterfort/cescel.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvis.h"
#include "asterfort/infniv.h"
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
    character(len=19) :: chanom
    character(len=*) :: ncmpva, ncmpvm
    character(len=8) :: nommod, nomaas, nomgd
    character(len=4) :: typech
    character(len=3) :: prolz
    character(len=8) :: crit, param
    character(len=24) :: option
    character(len=*) :: nochmd, nomamd
!
    integer(kind=8) :: nrofic, typent
    integer(kind=8) :: nbcmpv
    integer(kind=8) :: iinst, numpt, numord
    integer(kind=8) :: nbpgma(*), nbpgmm(*), nbspmm(*)
    integer(kind=8) :: codret
!
    real(kind=8) :: inst
    real(kind=8) :: prec
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='LRCEME')
!
    integer(kind=8) :: iaux, naux, unite, nbcham, nbcmp, jcmp
    integer(kind=8) :: tycha, junit
    integer(kind=8) :: vali(2)
    integer(kind=8) :: ibid
    integer(kind=8) :: jcesl
    med_idt :: idfimd
    integer(kind=8) :: ifm, nivinf, nseqca
    integer(kind=8) :: jnocmp, ncmprf, jcmpva, nbcmpa, iret, i, j, nncp
    integer(kind=8) :: edlect
    parameter(edlect=0)
!
    character(len=1) :: saux01
    character(len=8) :: saux08
    character(len=19) :: chames, ligrel
    character(len=64) :: nomcha
    character(len=64) :: valk(1)
    character(len=200) :: nofimd
    character(len=255) :: kfic
!
    aster_logical :: ttt
!
    call jemarq()
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'DEBUT DE '//nompro
    end if
1001 format(/, 10('='), a, 10('='),/)
!
!====
! 1. ALLOCATION D'UN CHAM_ELEM_S  (CHAMES)
!====
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
! 1.2. ==> ALLOCATION DU CHAM_ELEM_S
!
!               1234567890123456789
    chames = '&&      .CES.MED   '
    chames(3:8) = nompro
    ligrel = nommod//'.MODELE'
    if (nommod .eq. ' ') ligrel = ' '
!
    call jeexin(ncmpva, iret)
    if (iret .gt. 0) then
        call jeveuo(ncmpva, 'L', jcmpva)
        call jelira(ncmpva, 'LONMAX', nbcmpa)
        if (nomgd(1:4) .eq. 'VARI') then
            jnocmp = jcmpva
            ncmprf = nbcmpa
        else if (nbcmpa .le. ncmprf) then
            do i = 1, nbcmpa
                ttt = .false.
                do j = 1, ncmprf
                    if (zk8(jcmpva+i-1) .eq. zk8(jnocmp+j-1)) then
                        ttt = .true.
                    end if
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
!
        call getvis(' ', 'UNITE', scal=unite, nbret=iaux)
        call ulisog(unite, kfic, saux01)
        if (kfic(1:1) .eq. ' ') then
            call codent(unite, 'G', saux08)
            nofimd = 'fort.'//saux08
        else
            nofimd = kfic(1:200)
        end if
        call as_med_open(idfimd, nofimd, edlect, iret)
        call as_mfdnfd(idfimd, nbcham, iret)
        do i = 1, nbcham
            call as_mfdnfc(idfimd, i, nbcmp, iret)
            call wkvect('&&LRCEME.NOMCMP_K16', 'V V K16', nbcmp, jcmp)
            call wkvect('&&LRCEME.UNITCMP', 'V V K16', nbcmp, junit)
            call as_mfdfdi(idfimd, i, nomcha, tycha, zk16(jcmp), &
                           zk16(junit), nseqca, iret)
            if (nomcha .eq. nochmd) then
                ncmprf = nbcmp
                call wkvect('&&LRCEME.NOMCMP_K8', 'V V K8', nbcmp, jnocmp)
                do j = 1, nbcmp
                    zk8(jnocmp+j-1) = zk16(jcmp+j-1) (1:8)
                end do
                call jedetr('&&LRCEME.NOMCMP_K16')
                call jedetr('&&LRCEME.UNITCMP')
                call as_mficlo(idfimd, iret)
                goto 780
            end if
            call jedetr('&&LRCEME.NOMCMP_K16')
            call jedetr('&&LRCEME.UNITCMP')
        end do
        call as_mficlo(idfimd, iret)
    end if
!
780 continue
!
!====
! 3. LECTURE POUR CHAQUE TYPE DE SUPPORT
!====
!
    call lrcame(nrofic, nochmd, nomamd, nomaas, ligrel, &
                option, param, typech, typent, nbpgma, &
                nbpgmm, nbspmm, nbcmpv, ncmpva, ncmpvm, &
                iinst, numpt, numord, inst, crit, &
                prec, nomgd, ncmprf, jnocmp, chames, &
                codret)
!
    call jeveuo(chames//'.CESL', 'E', jcesl)
!
!====
! 4. TRANSFORMATION DU CHAM_ELEM_S EN CHAM_ELEM :
!====
!
!
    if (typech(1:4) .eq. 'CART') then
        call cescar(chames, chanom, 'V')
        nncp = 0
    else
        call cescel(chames, ligrel, option, param, prolz, &
                    nncp, 'V', chanom, 'F', ibid)
    end if
    if (nncp .gt. 0) then
        iaux = 0
        call jelira(chames//'.CESL', 'LONMAX', naux)
        do i = 1, naux
            if (zl(jcesl+i-1)) iaux = iaux+1
        end do
        vali(1) = iaux
        vali(2) = nncp
        valk(1) = nochmd
        call utmess('A', 'MED_83', nk=1, valk=valk, ni=2, &
                    vali=vali)
    end if
!
    call detrsd('CHAM_ELEM_S', chames)
!
!====
! 5. BILAN
!====
!
    if (codret .ne. 0) then
        call utmess('A', 'MED_55', sk=chanom)
    end if
!
!      IF(TYPECH(1:4).EQ.'ELGA')THEN
!        CALL JEDETR('&&LRCEME_NBPG_MAILLE')
!      ENDIF
    call jedetr('&&LRCEME.NOMCMP_K8')
    call jedema()
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'FIN DE '//nompro
    end if
!
end subroutine
