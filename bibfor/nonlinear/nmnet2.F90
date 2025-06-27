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

subroutine nmnet2(zimat, nmnbn, cnbn, cplas, czef, &
                  czeg, cief, cdeps, cdtg, cier, &
                  dc1, dc2, depsp2, normm)
    implicit none
!
!     CALCUL DE CNBN ET CDEPSP
!     QUAND DEUX CRITERES PLASTIQUES SONT ACTIVES
!
! IN  ZIMAT : ADRESSE DE LA LISTE DE MATERIAU CODE
! IN  NMNBN : FORCE - BACKFORCE
! IN  CDTG : MATRICE TANGENTE
! IN  DC1 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
! IN  DC2 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
!
! OUT CNBN : NOUVELLE FORCE - BACKFORCE
! OUT CPLAS : NOUVEAUX MOMENTS LIMITES DE PLASTICITE
! OUT CZEF : NOUVEAU ZERO ADIMENSIONNEL POUR LE CRITERE F
! OUT CZEG : NOUVEAU ZERO ADIMENSIONNEL POUR LE CRITERE G
! OUT CIEF : NOUVEAU CIEF > 0 : NBN HORS DE LA ZONE DE DEFINITION DE MP
! OUT CDEPS : NOUVEL INCREMENT DE DEFORMATION DANS LE REPERE ORTHO
! OUT CIER : NOUVEAU CODE ERREUR
! OUT DEPSP2 : NOUVEL INCREMENT DE DEF PLASTIQUE DANS LE REPERE ORTHO
! OUT NORMM : NORME SUR LA FONCTION MP = F(N)
!
#include "asterfort/gplass.h"
#include "asterfort/mppffn.h"
    integer(kind=8) :: j, cief, cier, zimat
!
    real(kind=8) :: nmnbn(6), cnbn(6), cplas(2, 3), czef, czeg
    real(kind=8) :: cdeps(6), cdtg(6, 6), dc1(6, 6), dc2(6, 6), normm
    real(kind=8) :: depsp2(6, 2), cp(6)
!
    cnbn = matmul(cdtg, cdeps)
    cp = matmul(dc1, depsp2(:, 1))
!
    do j = 1, 6
        cnbn(j) = nmnbn(j)+cnbn(j)-cp(j)
    end do
!
    cp = matmul(dc2, depsp2(:, 2))
!
    do j = 1, 6
        cnbn(j) = cnbn(j)-cp(j)
    end do
!
!     CALCUL DES MOMENTS LIMITES DE PLASTICITE
!     ET DES ZEROS DES CRITERES
    call mppffn(zimat, cnbn, cplas, czef, czeg, &
                cief, normm)
!
    if (cief .gt. 0) then
        cier = 2
        goto 30
    end if
!
    if ((gplass(cnbn, cplas, 1) .gt. czeg) .or. (gplass(cnbn, cplas, 2) .gt. czeg)) then
        cier = 1
    else
        cier = 0
    end if
!
30  continue
!
end subroutine
