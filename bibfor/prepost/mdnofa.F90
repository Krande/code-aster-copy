! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine mdnofa(numfam, nbfaex, nofaex, nomfam)
! person_in_charge: nicolas.sellenet at edf.fr
!  ENTREES :
!     NUMFAM = NUMERO DE LA FAMILLE
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
    integer :: numfam
    integer :: nbfaex
    character(len=*) :: nomfam
    character(len=*) :: nofaex(*)
!
! 0.3. ==> VARIABLES LOCALES
!
    integer :: jaux
    character(len=8) :: saux08
!
!
!====
! 1. ON INITIALISE AVEC DES ' '
!====
!
    nomfam = ' '
!
    nomfam(1:8) = 'FAMILLE_'
!
    call codent(numfam, 'G', saux08)
!
    jaux = lxlgut(saux08)
    nomfam(9:8+jaux) = saux08(1:jaux)
!
!====
! 2. MEMORISATION DU NOM DANS LE TABLEAU
!====
!
    nofaex(nbfaex+1) = nomfam
!
end subroutine
