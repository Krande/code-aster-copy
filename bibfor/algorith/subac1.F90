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
subroutine subac1(laxi, nno, vff, dff, geom, &
                  cova)
!
    implicit none
#include "asterf_types.h"
!
    aster_logical :: laxi
    integer(kind=8) :: nno
    real(kind=8) :: vff(nno), dff(nno), geom(2, nno), cova(3, 3)
!.......................................................................
!     CALCUL DE LA BASE COVARIANTE POUR UN ELEMENT LINEIQUE
!.......................................................................
! IN  AXI     TRUE SI AXI, FALSE SI D_PLAN OU C_PLAN
! IN  NNO     NOMBRE DE NOEUDS
! IN  VFF     VALEUR DES FONCTIONS DE FORMES
! IN  DFF     DERIVEE DES F. DE FORME
! IN  GEOM    COORDONNEES DES NOEUDS
! OUT COVA    COORDONNEES DES VECTEURS DE LA BASE COVARAINTE
!.......................................................................
!
    integer(kind=8) :: n, i
    real(kind=8) :: norme
!
    do i = 1, 3
        cova(i, 1) = 0.d0
        cova(i, 2) = 0.d0
    end do
!
!    CALCUL DU PREMIER VECTEUR TANGENT
    do n = 1, nno
        do i = 1, 2
            cova(i, 1) = cova(i, 1)+dff(n)*geom(i, n)
        end do
    end do
!
!    CALCUL DU SECOND VECTEUR TANGENT
    if (laxi) then
        do n = 1, nno
            cova(3, 2) = cova(3, 2)+vff(n)*geom(1, n)
        end do
    else
        cova(3, 2) = 1.d0
    end if
!
!    CALCUL DE LA NORMALE (PRODUIT VECTORIEL DES VECTEURS TANGENTS)
    cova(1, 3) = cova(2, 1)*cova(3, 2)-cova(3, 1)*cova(2, 2)
    cova(2, 3) = cova(3, 1)*cova(1, 2)-cova(1, 1)*cova(3, 2)
    cova(3, 3) = cova(1, 1)*cova(2, 2)-cova(2, 1)*cova(1, 2)
!
    norme = sqrt(cova(1, 3)**2+cova(2, 3)**2+cova(3, 3)**2)
    cova(1, 3) = cova(1, 3)/norme
    cova(2, 3) = cova(2, 3)/norme
    cova(3, 3) = cova(3, 3)/norme
!
end subroutine
