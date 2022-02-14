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

subroutine ut3mlg(nno, nnc, pgl, matl, matg)
    implicit none
!
#include "asterfort/assert.h"
!
    integer      :: nno, nnc
    real(kind=8) :: matl(*), pgl(3, 3), matg(*)
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
    integer         :: ii, jj, kk, nnoc, indx
    real(kind=8)    :: MatPassGL(nno*nnc, nno*nnc), mattmp(nno*nnc, nno*nnc)
! --------------------------------------------------------------------------------------------------
!
! Pour les index des matrices stockées en vecteur colonne
!   klv(i,j) : klv( indx(i,j,n) )
#define IndexMat(i,j,n) i+(j-1)*n
!
!   Pour l'instant cela ne marche qu'avec 3 composantes par noeuds
    ASSERT( nnc.eq.3 )
! --------------------------------------------------------------------------------------------------
    nnoc = nno*nnc
!   Matrice de passage Global vers Local
    MatPassGL(:,:) = 0.0
    if ( nno.eq.1 ) then
        do ii = 1, 3
            do jj = 1, 3
                MatPassGL(ii,jj) = pgl(ii,jj)
            enddo
        enddo
    else
        do ii = 1, 3
            do jj = 1, 3
                MatPassGL(ii  ,jj)   = pgl(ii,jj)
                MatPassGL(ii+3,jj+3) = pgl(ii,jj)
           enddo
        enddo
    endif
!
!   Matg = Pgl(^-1) * Matl * Pgl
!
!   Calcul de : mattmp() = Matl() * Pgl()
    mattmp(:,:) = 0.d0
    do ii = 1, nnoc
        do jj = 1, nnoc
            do kk = 1, nnoc
                indx = IndexMat(ii,kk,nnoc)
                mattmp(ii,jj) = mattmp(ii,jj) + matl( indx )*MatPassGL(kk,jj)
            enddo
        enddo
    enddo
!   Calcul de : matg() = transpose(Pgl) * mattmp()
    matg( nnoc*nnoc ) = 0.0
    do ii = 1, nnoc
        do jj = 1, nnoc
            indx = IndexMat(ii,jj,nnoc)
            do kk = 1, nnoc
                matg( indx ) = matg( indx ) + MatPassGL(kk,ii)*mattmp(kk,jj)
            enddo
        enddo
    enddo
!
end subroutine
