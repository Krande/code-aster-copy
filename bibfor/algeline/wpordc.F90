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
subroutine wpordc(type, shift, vp, x, m, &
                  neq)
    implicit none
#include "asterfort/utmess.h"
    integer(kind=8) :: type, neq, m
    complex(kind=8) :: x(neq, m), shift, vp(*)
!     TRI DES VALEURS (ET DES VECTEURS) PROPRES COMPLEXES
!     DEUX TYPE DE TRI :
!          - TRI DANS LE SPECTRE : SUIVANT ABS(SHIFT - VPQ)
!          - TRI DE PRESNTATION  : SUIVANT IM(VPQ) - IM(SHIFT)
!     ------------------------------------------------------------------
! IN  TYPE   : IS : TYPE DU TRI PAR ORDRE CROISSANT SUR LES VALEURS.
!                   * SI TYPE = 0  TRI DE PRESENTATION
!                   * SI TYPE = 1  TRI DANS LE SPECTRE
! IN  M      : IS : NOMBRE DE VALEUR PROPRE
! IN  SHIFT  : C8 : DECALAGE SPECTRAL
! VAR VP     : C8 : TABLEAU DES DES VALEURS PROPRES
! VAR X      : C8 : MATRICE DES VECTEURS PROPRES
! IN  NEQ    : IS : NOMBRE D'EQUATIONS
!                 SI NEQ < NBPRO ALORS ON NE TRIE PAS DE VECTEURS
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: p, om
    complex(kind=8) :: c, q
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    om = dimag(shift)
    if (type .eq. 0) then
        do i = 1, m, 1
            k = i
            p = dimag(vp(i))-om
            do j = i+1, m
                if ((dimag(vp(j))-om) .lt. p) then
                    p = dimag(vp(j))-om
                    k = j
                end if
            end do
            if (k .ne. i) then
                q = vp(i)
                vp(i) = vp(k)
                vp(k) = q
                do j = 1, neq, 1
                    c = x(j, i)
                    x(j, i) = x(j, k)
                    x(j, k) = c
                end do
            end if
        end do
    else if (type .eq. 1) then
        do i = 1, m, 1
            k = i
            p = abs(vp(i)-shift)
            do j = i+1, m
                if ((abs(vp(j)-shift)) .lt. p) then
                    p = abs(vp(j)-shift)
                    k = j
                end if
            end do
            if (k .ne. i) then
                q = vp(i)
                vp(i) = vp(k)
                vp(k) = q
                do j = 1, neq, 1
                    c = x(j, i)
                    x(j, i) = x(j, k)
                    x(j, k) = c
                end do
            end if
        end do
    else
        call utmess('F', 'ALGELINE3_97')
    end if
end subroutine
