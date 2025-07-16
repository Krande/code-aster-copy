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

subroutine ut3mlg(nno, nnc, pgl, matl, matg)
!
    use linalg_ops_module, only: as_matmul
!
    implicit none
!
#include "asterfort/assert.h"
!
    integer(kind=8)      :: nno, nnc
    real(kind=8) :: matl(nno*nnc, nno*nnc), pgl(3, 3), matg(nno*nnc, nno*nnc)
!
! --------------------------------------------------------------------------------------------------
!
!           TRANSFORMATION DES MATRICES ELEMENTAIRES NON-SYMÉTRIQUES
!                   PASSAGE DU REPERE LOCAL AU REPERE GLOBAL
!
!   In
!       nno     nombre de noeuds
!       nnc     nombre de composantes
!       pgl     matrice de passage local -> global
!       matl    matrice complète dans le repère local (stockage en colonne)
!
!   Out
!       matg    matrice dans le repère global (stockage en colonne)
!
! --------------------------------------------------------------------------------------------------
!
!   !!!! Les matrices sont complètes et stockées en colonnes !!!!
!           nombre d'élément des matrices : (nno*nnc)*(nno*nnc)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)         :: ii, jj
    real(kind=8)    :: MatPassGL(nno*nnc, nno*nnc)
! --------------------------------------------------------------------------------------------------
!
!   Pour l'instant cela ne marche qu'avec 3 composantes par noeuds
    ASSERT(nnc .eq. 3)
! --------------------------------------------------------------------------------------------------
!   Matrice de passage Global vers Local
    MatPassGL(:, :) = 0.0
    if (nno .eq. 1) then
        do ii = 1, 3
            do jj = 1, 3
                MatPassGL(ii, jj) = pgl(ii, jj)
            end do
        end do
    else
        do ii = 1, 3
            do jj = 1, 3
                MatPassGL(ii, jj) = pgl(ii, jj)
                MatPassGL(ii+3, jj+3) = pgl(ii, jj)
            end do
        end do
    end if

    matg = as_matmul(transpose(MatPassGL), as_matmul(Matl, MatPassGL))

end subroutine
