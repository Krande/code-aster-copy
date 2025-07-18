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
subroutine gdfint(kp, nno, ajacob, pjacob, en, &
                  enprim, x0pg, pn, pm, fint)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           CONTRIBUTION DU POINT DE GAUSS NUMERO KP AUX FORCES INTERNES
!
!     IN  : KP        : NUMERO DU POINT DE GAUSS
!           NNO       : NOMBRE DE NOEUDS DE L'ELEMENT
!           AJACOB    : JACOBIEN
!           PJACOB    : POIDS * JACOBIEN
!           EN        : FONCTIONS DE FORME
!           ENPRIM    : DERIVEES DES FONCTIONS DE FORME
!           X0PG      : DERIVEES DES COORDONNEES PAR RAP. A L'ABS. CURV.
!           PN        : RESULTANTE DES FORCES AU PT DE GAUSS EN AX.GENE.
!           PM        : MOMENT RESULTANT AU PT DE GAUSS EN AXES GENERAUX
!
!     OUT : FINT      : FORCES INT. (CUMUL DES CONTRIB. DES PTS DE GAUS)
! ------------------------------------------------------------------
    implicit none
#include "asterfort/gdmb.h"
#include "asterfort/promat.h"
#include "asterfort/transp.h"
    real(kind=8) :: en(3, 2), enprim(3, 2), x0pg(3), pn(3), pm(3), fint(6, 3)
    real(kind=8) :: b(6, 6), bt(6, 6), vect(6), fors(6)
    integer(kind=8) :: i, k, kp, ne, nno
    real(kind=8) :: ajacob, pjacob
!-----------------------------------------------------------------------
!
    do ne = 1, nno
        call gdmb(ne, kp, ajacob, en, enprim, &
                  x0pg, b)
        call transp(b, 6, 6, 6, bt, &
                    6)
        do i = 1, 3
            vect(i) = pn(i)
            vect(3+i) = pm(i)
        end do
        call promat(bt, 6, 6, 6, vect, &
                    6, 6, 1, fors)
        do k = 1, 6
            fint(k, ne) = fint(k, ne)+pjacob*fors(k)
        end do
    end do
end subroutine
