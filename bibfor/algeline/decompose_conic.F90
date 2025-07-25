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

subroutine decompose_conic(m0, nline, line1, line2, indic)
!
    implicit none
#include "asterfort/mat_com.h"
#include "asterfort/num_rank_mat33.h"
!
    real(kind=8), intent(in) :: m0(3, 3)
    integer(kind=8), intent(out) :: nline
    real(kind=8), intent(out) :: line1(3)
    real(kind=8), intent(out) :: line2(3)
    real(kind=8), intent(out), optional :: indic
!
!
!    TROUVER LES DROITES CORRESPONDANTS A UNE CONIQUE DEGENEREE
!    EN UTILISANT LA GEOMETRIE PROJECTIVE (COORDONNES HOMOGENES)
!
! IN   M0  : MATRICE TAILLE 3*3 DE LA CONIQUE DEGENEREE
! OUT  NLINE : NOMBRE DE LIGNES TROUVEES
! OUT  LINE1 : COORDONNEES HOMOGENES DE LA PREMIERE LIGNE TROUVEE
! OUT  LINE2 : COORDONNEES HOMOGENES DE LA DEUXIEME LIGNE TROUVEE
! OUT  INDIC : DISTANCE AU CENTRE DU POINT D INTERSECTION ENTRE LES LIGNES
!  PLUS IL EST PROCHE DU CENTRE, MIEUX C'EST POUR LA RESOLUTION
!  (ON VEUT EVITER LES DROITES QUASIMENT PARALLELES)
!
    real(kind=8) :: prec, mat_b(3, 3), matrice(3, 3)
    real(kind=8) :: max_di, b, p(3), Mp(3, 3), test, val_max, condit
    integer(kind=8) :: rank_m0, imax, jmax, i, j
!
!     Initialisations
    line1 = (/0.d0, 0.d0, 0.d0/)
    line2 = (/0.d0, 0.d0, 0.d0/)
    prec = 1.e-13
    rank_m0 = num_rank_mat33(m0, prec, condit)
    if (rank_m0 .lt. 2) then
        matrice = m0
        ! deux droites paralleles
        ! on prend en dernier recours seulement
        indic = 1.d20
    else
        mat_b = -mat_com(3, m0)
        !
        max_di = 0.d0
        imax = 0
        !
        ! mat_b(3,3) doit etre > 0 pour lignes s intersectant
        ! en fait tous les coefficients
        do i = 1, 3
            if (abs(mat_b(i, i)) .gt. max_di) then
                max_di = abs(mat_b(i, i))
                imax = i
            end if
        end do
        !
        ! Verification
        ! Comme on a pris le max, on peut se contenter de regarder par rapport a 0
        !
        ! on peut ajouter un test sur b(3,3) pour etre sur d avoir des droites secantes
        if (mat_b(imax, imax) .le. 0.d0 .or. &
            mat_b(3, 3) .le. 0.d0) then
            nline = 0
            goto 99
        end if
        !
        ! on peut ajouter un test sur b(3,3) pour etre sur d avoir des droites secantes
        b = sqrt(mat_b(imax, imax))
        p(:) = mat_b(:, imax)/b
        !
        ! matrice produit vectoriel
        do i = 1, 3
            Mp(i, i) = 0.d0
        end do
        Mp(1, 2) = p(3)
        Mp(1, 3) = -p(2)
        Mp(2, 1) = -p(3)
        Mp(2, 3) = p(1)
        Mp(3, 1) = p(2)
        Mp(3, 2) = -p(1)
        !
        !
        matrice = m0+Mp
        !
        indic = sqrt((p(1)*p(1)+p(2)*p(2))/(p(3)*p(3)))
    end if
!
!     La matrice etant de rang 1
!     On recupere la valeur max de toute la matrice
    val_max = 0.d0
    imax = 0
    jmax = 0
    do i = 1, 3
        do j = 1, 3
            test = abs(matrice(i, j))
            if (test .gt. val_max) then
                val_max = test
                imax = i
                jmax = j
            end if
        end do
    end do
!
!     On construit la decomposition ligne colonne a partir de cet element max
    nline = 2
    line1(:) = matrice(imax, :)
    line2(:) = matrice(:, jmax)
99  continue
!
end subroutine
