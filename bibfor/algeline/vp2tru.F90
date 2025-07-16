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

subroutine vp2tru(method, ty, alpha, beta, signes, &
                  a, nbvect, w, z, wk, &
                  mxiter, ier, nitqr)
    implicit none
#include "asterfort/utmess.h"
#include "asterfort/vphqrp.h"
    character(len=1) :: ty
    character(len=8) :: method
    integer(kind=8) :: nbvect, mxiter, ier, nitqr
    integer(kind=8) :: vali
    real(kind=8) :: alpha(nbvect), beta(nbvect), signes(nbvect)
    real(kind=8) :: a(nbvect, *), w(*), z(*), wk(*)
    real(kind=8) :: valr(2)
!     EXTENSION DE LA METHODE DE LANCZOS POUR TRIDIAGONALISER UN
!     FAISCEAU DE MATRICE INDEFINI
!     ------------------------------------------------------------------
! IN  METHOD : K8 : METHODE DE RESOLUTION
!         SI 'TRI_DIAG' : A EST CONSTRUITE A PARTIR DE ALPHA,BETA,SIGNES
!         SI 'ARNOLDI'  : A EST DEJA CONSTRUITE
! IN  TY : K1 : TYPE DE PROBLEME TRAITE A L'ORIGINE
!               'G' GENERALISE  ==> ELEMENTS PROPRES REELS
!               'Q' QUADRATIQUE ==> ELEMENTS PROPRES COMPLEXES
! IN  ALPHA : DIAGONALE DE LA MATRICE TRIDIAGONALE ('TRI_DIAG')
! IN  BETA  : SURDIAGONALE DE LA MATRICE TRIDIAGONALE ('TRI_DIAG')
! IN SIGNES:SIGNE POUR PASSER DE LA SUR A LA SOUS DIAGONALE('TRI_DIAG')
! OUT W : C : VALEURS PROPRES DU SYSTEME
! OUT ALPHA : R : PARTIE REELLE DES VALEURS PROPRES DU SYSTEME
! OUT BETA  : R : PARTIE IMAGINAIRE DES VALEURS PROPRES DU SYSTEME
! OUT Z : C : VECTEURS PROPRES DU SYSTEME
! OUT NITQR : NOMBRE D'ITERATIONS QR POUT ATTEIENDRE LA CONVERGENCE
!
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ival1, j
    real(kind=8) :: ww
!-----------------------------------------------------------------------
    if (method .eq. 'TRI_DIAG') then
        a(:, 1:nbvect) = 0.d0
        a(1, 1) = signes(1)*alpha(1)
        a(1, 2) = signes(1)*beta(2)
        do i = 2, nbvect-1
            a(i, i-1) = signes(i)*beta(i)
            a(i, i) = signes(i)*alpha(i)
            a(i, i+1) = signes(i)*beta(i+1)
        end do
        a(nbvect, nbvect-1) = signes(nbvect)*beta(nbvect)
        a(nbvect, nbvect) = signes(nbvect)*alpha(nbvect)
    end if
!
    ival1 = 1
    ier = 0
    call vphqrp(a, nbvect, nbvect, ival1, w, &
                z, nbvect, wk, mxiter, ier, &
                nitqr)
    if (ier .ne. 0) then
        call utmess('F', 'ALGELINE3_57')
    end if
!
    if (ty .eq. 'G') then
        do i = 1, nbvect
            alpha(i) = w(2*i-1)
            if (abs(w(2*i)) .gt. 1.d-75) then
                ww = w(2*i)/w(2*i-1)
                vali = i
                valr(1) = w(2*i-1)
                valr(2) = ww
                call utmess('I', 'ALGELINE4_65', si=vali, nr=2, valr=valr)
            end if
            beta(i) = 0.d0
            do j = 1, nbvect
                a(i, j) = z(2*nbvect*(j-1)+2*i-1)
            end do
        end do
    else
        do i = 1, nbvect
            alpha(i) = w(2*i)
            beta(i) = w(2*i-1)
        end do
    end if
!
end subroutine
