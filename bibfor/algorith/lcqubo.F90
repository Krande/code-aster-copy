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

subroutine lcqubo(ep0, ep1, l0, l1, etamin, &
                  etamax, vide, etam, etap)
    implicit none
#include "asterf_types.h"
    real(kind=8), intent(in) :: ep0(6), ep1(6), l0, l1, etamin, etamax
    real(kind=8), intent(out) :: etam, etap
    aster_logical, intent(out) :: vide
!
! ----------------------------------------------------------------------
!  BORNES POUR LE PILOTAGE RELATIF AU CRITERE QUADRATIQUE
!  ON NE FAIT RIEN A L'HEURE ACTUELLE
!    CHI**2 + L0 + L1*ETA = 0
! ----------------------------------------------------------------------
!  IN  EP0    DEFORMATION FIXE
!  IN  EP1    DEFORMATION PILOTEE
!  IN  L0,L1  COMPOSANTES DU TERME AFFINE
!  IN  ETAMIN BORNE MIN INITIALE
!  IN  ETAMAX BORNE MAX INITIALE
!  OUT VIDE   CODE RETOUR: T->PAS DE SOLUTION ; F->BORNES
!  OUT ETAM   NOUVELLE BORNE MIN
!  OUT ETAP   NOUVELLE BORNE MAX
! ----------------------------------------------------------------------
    etam = etamin
    etap = etamax
    vide = etamax .le. etamin
end subroutine
