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

subroutine rssepa(resultZ, numeStore, modelZ, materFieldZ, caraElemZ, listLoadZ)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/rsadpa.h"

    character(len=*), intent(in) :: resultZ
    integer(kind=8), intent(in)  :: numeStore
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, listLoadZ
!----------------------------------------------------------------------
!     BUT: ECRIRE DANS LA SD RESULTAT LES PARAMETRES MODELE, MATE,
!          CARELE ET EXCIT POUR LE NUME_ORDRE NUORDR
!
!     IN      RESULT : NOM DE LA SD RESULTAT
!     IN      IORDR  : NUMERO D'ORDRE
!     IN      MODELE : NOM DU MODELE
!     IN      MATE   : NOM DU CHAMP MATERIAU
!     IN      CARELE : NOM DE LA CARACTERISTIQUE ELEMENTAIRE
!     IN      EXCIT  : NOM DE LA SD INFO_CHARGE
!
!
    character(len=19) :: listLoad
    integer(kind=8) :: jvPara
! ----------------------------------------------------------------------
!
    listLoad = listLoadZ(1:19)
    call rsadpa(resultZ, 'E', 1, 'MODELE', numeStore, 0, sjv=jvPara)
    zk8(jvPara) = modelZ(1:8)
    call rsadpa(resultZ, 'E', 1, 'CHAMPMAT', numeStore, 0, sjv=jvPara)
    zk8(jvPara) = materFieldZ(1:8)
    call rsadpa(resultZ, 'E', 1, 'CARAELEM', numeStore, 0, sjv=jvPara)
    zk8(jvPara) = caraElemZ(1:8)
    call rsadpa(resultZ, 'E', 1, 'EXCIT', numeStore, 0, sjv=jvPara)
    zk24(jvPara) = listLoad
!
end subroutine
