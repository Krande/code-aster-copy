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
subroutine mdnofa(numfam, nogrf, nbgf, nbfaex, nofaex, &
                  nomfam)
! person_in_charge: nicolas.sellenet at edf.fr
!  ENTREES :
!     NUMFAM = NUMERO DE LA FAMILLE
!     NOGRF  = NOMS DES GROUPES D ENTITES DE LA FAMILLE
!     NBGF   = NOMBRE DE GROUPES DE LA FAMILLE
!     NBFAEX = NOMBRE DE FAMILLES EXISTANTES POUR CETTE SERIE
!     NOFAEX = NOMS DES FAMILLES DEJA CREEES
!  SORTIES :
!     NOMFAM = NOM DE LA FAMILLE
!     ------------------------------------------------------------------
!
!    LE NOM DE LA FAMILLE SERA LA CONCATENATION DES NOMS DES GROUPES.
!    ON PLACE UN '_' ENTRE DEUX NOMS
!    EXEMPLE :
!       GROUPE 1 : 'HAUT'
!       GROUPE 2 : 'FACE_X'
!       FAMILLE ==> 'HAUT_FACE_X'
!                    12345678901234567890123456789012
!
!     SI ON N'A PAS LA PLACE, ON MET UN NOM ARBITRAIRE
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterfort/codent.h"
#include "asterfort/lxlgut.h"
    integer(kind=8) :: numfam
    integer(kind=8) :: nbgf
    integer(kind=8) :: nbfaex
    character(len=*) :: nogrf(*)
    character(len=*) :: nomfam
    character(len=*) :: nofaex(*)
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    integer(kind=8) :: iaux, jaux
    integer(kind=8) :: ideb, ifin
    integer(kind=8) :: lgnofa, lgnofx
!
    character(len=8) :: saux08
!
!
!
!====
! 1. ON INITIALISE AVEC DES ' '
!====
!
    lgnofx = len(nomfam)
    do iaux = 1, lgnofx
        nomfam(iaux:iaux) = ' '
    end do
!
!====
! 2. CALCUL DE LA TAILLE NECESSAIRE A LA CONCATENATION
!====
!
    lgnofa = 0
    do iaux = 1, nbgf
        lgnofa = lgnofa+lxlgut(nogrf(iaux))
        if (iaux .lt. nbgf) then
            lgnofa = lgnofa+1
        end if
    end do
!
!====
! 2. ON A LA PLACE : FABRICATION DU NOM
!====
!
    if (lgnofa .le. lgnofx) then
!
        ifin = 0
        do iaux = 1, nbgf
            jaux = lxlgut(nogrf(iaux))
            ideb = ifin+1
            ifin = ifin+jaux
            nomfam(ideb:ifin) = nogrf(iaux) (1:jaux)
            if (iaux .lt. nbgf) then
                ifin = ifin+1
                nomfam(ifin:ifin) = '_'
            end if
        end do
!
!====
! 4. SINON, C'EST UN NOM ARBITRAIRE, CONSTRUIT AVEC LE NUMERO DE
!    LA FAMILLE
!====
!
    else
!
!                      12345678
        nomfam(1:8) = 'FAMILLE_'
!
        call codent(numfam, 'G', saux08)
!
        jaux = lxlgut(saux08)
        nomfam(9:8+jaux) = saux08(1:jaux)
!
    end if
!
!====
! 5. MEMORISATION DU NOM DANS LE TABLEAU
!====
!
    nofaex(nbfaex+1) = nomfam
!
end subroutine
