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

subroutine ircnme(ifi, nochmd, chanom, typech, modele, &
                  nbcmp, nomcmp, partie, numpt, instan, &
                  numord, nbnoec, linoec, sdcarm, carael, &
                  field_type, lfichUniq, codret)
!_______________________________________________________________________
! person_in_charge: nicolas.sellenet at edf.fr
!        IMPRESSION DU CHAMP CHANOM NOEUD ENTIER/REEL
!        AU FORMAT MED
!     ENTREES:
!       IFI    : UNITE LOGIQUE D'IMPRESSION DU CHAMP
!       NOCHMD : NOM MED DU CHAM A ECRIRE
!       PARTIE: IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!               UN CHAMP COMPLEXE
!       CHANOM : NOM ASTER DU CHAM A ECRIRE
!       TYPECH : TYPE DU CHAMP
!       MODELE : MODELE ASSOCIE AU CHAMP
!       NBCMP  : NOMBRE DE COMPOSANTES A ECRIRE
!       NOMCMP : NOMS DES COMPOSANTES A ECRIRE
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       NBNOEC : NOMBRE DE NOEUDS A ECRIRE (O, SI TOUS LES NOEUDS)
!       LINOEC : LISTE DES NOEUDS A ECRIRE SI EXTRAIT
!       SDCARM : CARA_ELEM (UTILE POUR LES SOUS-POINTS)
! In  field_type       : type of field (symbolic name in result datastructure)
!    SORTIES:
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_______________________________________________________________________
!
    implicit none
#include "jeveux.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/ircame.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
! 0.1. ==> ARGUMENTS
!
    character(len=8) :: typech, modele, sdcarm, carael
    character(len=19) :: chanom
    character(len=64) :: nochmd
    character(len=*) :: nomcmp(*), partie
!
    integer(kind=8) :: nbcmp, numpt, ifi, numord
    integer(kind=8) :: nbnoec
    integer(kind=8) :: linoec(*)
!
    real(kind=8) :: instan
    character(len=16), intent(in) :: field_type
!
    aster_logical :: lfichUniq
!
    integer(kind=8) :: codret
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='IRCNME')
!
    character(len=19) :: chamns
!
    integer(kind=8) :: jcnsk, jcnsd, jcnsc, jcnsv, jcnsl, nbCmpDyna
!     ------------------------------------------------------------------
!
    call jemarq()
!
!====
! 1. PREALABLE
!====
!
!    --- CONVERSION CHAM_NO -> CHAM_NO_S
!               1234567890123456789
    chamns = '&&      .CNS.MED   '
    chamns(3:8) = nompro
    call cnocns(chanom, 'V', chamns)
!
!    --- ON RECUPERE LES OBJETS
!
    call jeveuo(chamns//'.CNSK', 'L', jcnsk)
    call jeveuo(chamns//'.CNSD', 'L', jcnsd)
    call jeveuo(chamns//'.CNSC', 'L', jcnsc)
    call jeveuo(chamns//'.CNSV', 'L', jcnsv)
    call jeveuo(chamns//'.CNSL', 'L', jcnsl)
!
!====
! 2. ECRITURE DES CHAMPS AU FORMAT MED
!====
!
    call ircame(ifi, nochmd, chanom, typech, modele, &
                nbcmp, nomcmp, ' ', partie, numpt, &
                instan, numord, jcnsk, jcnsd, jcnsc, &
                jcnsv, jcnsl, nbnoec, linoec, sdcarm, &
                carael, field_type, nbCmpDyna, lfichUniq, codret)
!
!====
! 3. ON NETTOIE
!====
!
    call detrsd('CHAM_NO_S', chamns)
!
!====
! 4. BILAN
!====
!
    if (codret .ne. 0 .and. codret .ne. 100) then
        call utmess('A', 'MED_89', sk=chanom)
    end if
!
    call jedema()
!
end subroutine
