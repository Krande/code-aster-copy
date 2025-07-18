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

subroutine mdexch(nofimd, idfimd, nochmd, numpt, numord, &
                  nbcmpc, nomcmc, nbvato, typent, typgeo, &
                  existc, nbcmfi, nmcmfi, nbval, nbprof, &
                  codret)
! person_in_charge: nicolas.sellenet at edf.fr
!_____________________________________________________________________
!        FORMAT MED : EXISTENCE D'UN CHAMP DANS UN FICHIER
!               - -   --             --
!_______________________________________________________________________
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       IDFIMD : OU NUMERO DU FCHIER MED DEJA OUVERT
!       NOCHMD : NOM MED DU CHAMP A CONTROLER
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       NUMORD : NUMERO D'ORDRE
!       NBCMPC : NOMBRE DE COMPOSANTES A CONTROLER
!       NOMCMC : SD DES NOMS DES COMPOSANTES A CONTROLER (K16)
!       NBVATO : NOMBRE DE VALEURS TOTAL
!       TYPENT : TYPE D'ENTITE MED DU CHAMP A CONTROLER
!       TYPGEO : TYPE GEOMETRIQUE MED DU CHAMP A CONTROLER
!     SORTIES:
!       EXISTC : 0 : LE CHAMP EST INCONNU DANS LE FICHIER
!               >0 : LE CHAMP EST CREE AVEC :
!                1 : LES COMPOSANTES VOULUES NE SONT PAS TOUTES
!                    ENREGISTREES
!                2 : AUCUNE VALEUR POUR CE TYPE ET CE NUMERO D'ORDRE
!                3 : DES VALEURS A CE NUMERO D'ORDRE
!                4 : DES VALEURS A CE NUMERO D'ORDRE, MAIS EN NOMBRE
!                    DIFFERENT
!       NBCMFI : NOMBRE DE COMPOSANTES DANS LE FICHIER
!       NMCMFI : SD DU NOM DES COMPOSANTES DANS LE FICHIER
!       NBVAL  : NOMBRE DE VALEURS DANS LE FICHIER
!       NBPROF : NOMBRE DE PROFILS PRESENTS
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_______________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/mdexcc.h"
#include "asterfort/mdexcv.h"
#include "asterfort/utmess.h"
    character(len=*) :: nofimd
    character(len=*) :: nochmd
    character(len=*) :: nomcmc, nmcmfi
!
    integer(kind=8) :: numpt, numord, nbcmpc, nbprof
    med_idt :: idfimd
    integer(kind=8) :: nbvato, typent, typgeo
    integer(kind=8) :: existc, nbcmfi, nbval
!
    integer(kind=8) :: codret
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
!
    integer(kind=8) :: iaux
!
!====
! 1. LE CHAMP A-T-IL ETE CREE ?
!====
!
    call mdexcc(nofimd, idfimd, nochmd, nbcmpc, nomcmc, &
                iaux, nbcmfi, nmcmfi, codret)
!
!====
! 2. SUITE DU DIAGNOSTIC
!====
!
! 2.1. ==> LE CHAMP N'EST PAS CREE
!
    if (iaux .eq. 0) then
!
        existc = 0
!
! 2.2. ==> LES COMPOSANTES VOULUES NE SONT PAS TOUTES ENREGISTREES
!
    else if (iaux .eq. 2) then
!
        existc = 1
!
! 2.3. ==> SI LE CHAMP EST CORRECTEMENT CREE, COMBIEN A-T-IL DE VALEURS
!          A CE COUPLE NUMERO DE PAS DE TEMPS / NUMERO D'ORDRE ?
!
    else if (iaux .eq. 1) then
!
        call mdexcv(nofimd, idfimd, nochmd, numpt, numord, &
                    typent, typgeo, nbval, nbprof, codret)
!
        if (nbval .eq. 0) then
            existc = 2
        else if (nbval .eq. nbvato) then
            existc = 3
        else
            existc = 4
        end if
!
! 2.4. ==> BIZARRE
!
    else
!
        call utmess('F', 'MED_76')
!
    end if
!
end subroutine
