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
subroutine utpnlg(nno, nnc, pgl, matl, mate)
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/r8inir.h"
!
    integer(kind=8) :: nno, nnc
    real(kind=8) :: mate(1), pgl(3, 3), matl(nno*nnc, nno*nnc)
! .....................................................................C
! .....................................................................C
!    TRANSFORMATION DES MATRICES ELEMENTAIRES                          C
!    PASSAGE DU REPERE LOCAL AU REPERE GLOBAL                          C
!                                                                      C
!    - ARGUMENTS:                                                      C
!        DONNEES:      NNO     -->  NOMBRE DE NOEUDS                   C
!                      NNC     -->  NOMBRE DE COMPOSANTES              C
!                      PGL     -->  MATRICE DE PASSAGE GLOBAL -> LOCAL C
!                      MATL    -->  MATRICE ELEMENTAIRE LOCALE         C
!                                                                      C
!        SORTIE :      MATE    -->  MATRICE ELEMENTAIRE GLOBALE        C
! .....................................................................C
! .....................................................................
    integer(kind=8) :: i, j, k, ii, nj
    real(kind=8) :: mt(nno*nnc, nno*nnc), matg(nno*nnc, nno*nnc)
! .....................................................................
    ASSERT(nnc .eq. 3)
    nj = nno*nnc
! --- MATRICE DE TRANSFERT
    call r8inir(nno*nno*nnc*nnc, 0.d0, mt, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 0, nno-1
                mt(i+k*3, j+k*3) = pgl(i, j)
            end do
        end do
    end do
!
! --- ON EFFECTUE : MATG() = MATL() * MT()
    do k = 1, nno*nnc
        do i = 1, nno*nnc
            matg(i, k) = 0.d0
            do j = 1, nj
                matg(i, k) = matg(i, k)+matl(i, j)*mt(j, k)
            end do
        end do
    end do
! --- MULTIPLICATION PAR LA MATRICE TRANSPOSEE DE "MT" LORSQUE
!           "MATE" EST RECTANGULAIRE DE DIMENSIONS 7X7
    do i = 1, nno*nnc
        ii = nj*(i-1)
        do k = 1, nno*nnc
            mate(ii+k) = 0.d0
            do j = 1, nj
                mate(ii+k) = mate(ii+k)+mt(j, i)*matg(j, k)
            end do
        end do
    end do
! .....................................................................
end subroutine
