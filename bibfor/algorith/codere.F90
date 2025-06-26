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
subroutine codere(cod, npg, codret)

    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"

    integer(kind=8) :: npg, cod(npg), codret
! --------------------------------------------------------------------------------------------------
!     SYNTHESE DES CODES RETOURS : EN ENTREE, ON A UN TABLEAU
!     DE DIM. NPG CONTENANT LES CODES RETOURS DE TOUS LES PTS DE
!     GAUSS. EN SORTIE, ON A UN SEUL CODE RETOUR RESUME
! --------------------------------------------------------------------------------------------------
!     in  cod     : tableau contenant les codes retours de tous les pts de gauss
!     in  npg     : nbre de pts de gauss de l'elelemt traite
!     out codret  : code retour resume
!         codret=0 : tout est ok
!         codret=1 : echec integration loi de comportement
!         codret=2 : resultats non valides physiquement
!         codret=3 : c_plan deborst sigzz non nul
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: NBR_CODRET = 4
    integer(kind=8), parameter, dimension(0:NBR_CODRET):: codret_to_gravite = [0, 4, 2, 3, 1]
    integer(kind=8), parameter, dimension(0:NBR_CODRET):: gravite_to_codret = [0, 4, 2, 3, 1]
    integer(kind=8) :: i, grave
! --------------------------------------------------------------------------------------------------

    ! Controle de l'etendue des codes retour
    ASSERT(maxval(cod) .le. NBR_CODRET .and. minval(cod) .ge. 0)

    ! Gravite maximale des codes retour : la gravite la plus elevee l'emporte sur les autres
    grave = 0
    do i = 1, npg
        grave = max(grave, codret_to_gravite(cod(i)))
    end do

    ! Code emis le plus grave
    codret = gravite_to_codret(grave)

end subroutine
