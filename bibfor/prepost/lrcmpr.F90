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

subroutine lrcmpr(idfimd, nomprf, ntproa, lgproa, codret)
!_____________________________________________________________________
! person_in_charge: nicolas.sellenet at edf.fr
! ======================================================================
!     LECTURE D'UN CHAMP - FORMAT MED - PROFIL
!     -    -       -              -     --
!-----------------------------------------------------------------------
!      ENTREES:
!       IDFIMD : IDENTIFIANT DU FICHIER MED
!       NOMPRF : NOM MED DU PROFIL A LIRE
!      SORTIES:
!       NTPROA : TABLEAU QUI CONTIENT LE PROFIL ASTER
!       LGPROA : LONGUEUR DU PROFIL ASTER
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_____________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/as_mpfprr.h"
#include "asterfort/as_mpfpsn.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    med_idt :: idfimd
    integer(kind=8) :: lgproa
    integer(kind=8) :: codret
!
    character(len=*) :: nomprf
    character(len=*) :: ntproa
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    character(len=8) :: saux08
    parameter(nompro='LRCMPR')
!
    integer(kind=8) :: ifm, nivinf
!
    integer(kind=8) :: adproa, adprom
    integer(kind=8) :: lgprom
    integer(kind=8) :: iaux
!
    character(len=24) :: ntprom
!
!====
! 1. PREALABLES
!====
!
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'DEBUT DE '//nompro
    end if
1001 format(/, 10('='), a, 10('='),/)
!
! 1.2. ==> NOMS DES TABLEAUX DE TRAVAIL
!               12   345678   9012345678901234
    ntprom = '&&'//nompro//'.PROFIL_MED     '
!
!====
! 2. NOMBRE DE VALEURS LIEES AU PROFIL
!====
!
    call as_mpfpsn(idfimd, nomprf, lgprom, codret)
    if (codret .ne. 0) then
        saux08 = 'mpfpsn'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    if (nivinf .gt. 1) then
        write (ifm, 4101) nomprf, lgprom
    end if
4101 format('. LECTURE DU PROFIL : ', a,&
   &     /, '... LONGUEUR : ', i8)
!
!====
! 3. LECTURE DES VALEURS DU PROFIL MED
!====
!
    call wkvect(ntprom, 'V V I', lgprom, adprom)
!
    call as_mpfprr(idfimd, zi(adprom), lgprom, nomprf, codret)
    if (codret .ne. 0) then
        saux08 = 'mpfprr'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    if (nivinf .gt. 1) then
        if (lgprom .ge. 10) then
            write (ifm, 4201) zi(adprom), zi(adprom+1), zi(adprom+2)
            write (ifm, 4202) zi(adprom+lgprom-3), zi(adprom+lgprom-2), &
                zi(adprom+lgprom-1)
        else
            write (ifm, 4203) (zi(adprom+iaux), iaux=0, lgprom-1)
        end if
    end if
4201 format('... 3 1ERES VALEURS     : ', 3i8)
4202 format('... 3 DERNIERES VALEURS : ', 3i8)
4203 format('... VALEURS : ', 10i8)
!
!====
! 4. TRANSFERT EN UN PROFIL ASTER
!====
!          EN FAIT, DANS LE CAS DES NOEUDS, IL Y A IDENTITE ENTRE LES
!          DEUX CAR ON NE RENUMEROTE PAS LES NOEUDS (CF IRMMNO)
!
    lgproa = lgprom
    call wkvect(ntproa, 'V V I', lgproa, adproa)
!
    do iaux = 0, lgprom-1
        zi(adproa+iaux) = zi(adprom+iaux)
    end do
!
    call jedetr(ntprom)
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'FIN DE '//nompro
    end if
!
end subroutine
