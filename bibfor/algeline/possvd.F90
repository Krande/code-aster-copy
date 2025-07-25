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
subroutine possvd(nm, m, n, w, matu, &
                  u, matv, v, eps, rg, &
                  rv1)
    implicit none
!
! DESCRIPTION :   POST-TRAITEMENTS AU CALCUL DE LA DECOMPOSITION AUX
! -----------     VALEURS SINGULIERES
!                                      T
!                             A = U S V
!
!                 REALISE PAR LA ROUTINE CALSVD
!                 A EST UNE MATRICE REELLE RECTANGULAIRE (M,N)
!
! IN     : NM   : INTEGER , SCALAIRE
!                 PREMIERE DIMENSION DES TABLEAUX A, U ET V, DECLAREE
!                 DANS L'APPELANT, NM >= MAX(M,N)
! IN     : M    : INTEGER , SCALAIRE
!                 NOMBRE DE LIGNES DES MATRICES A ET U
! IN     : N    : INTEGER , SCALAIRE
!                 NOMBRE DE COLONNES DES MATRICES A ET U
!                  = ORDRE DE LA MATRICE V
! IN/OUT : W    : REAL*8 , VECTEUR DE DIMENSION N
!                 CONTIENT LES N VALEURS SINGULIERES DE A
!                 EN SORTIE LES VALEURS SINGULIERES SONT REORDONNEES
!                 PAR MODULE DECROISSANT
! IN     : MATU : LOGICAL , SCALAIRE
!                 SI MATU = .TRUE. LA MATRICE U A ETE CALCULEE
! IN/OUT : U    : REAL*8 , TABLEAU DE DIMENSION (NM,N)
!                 SI MATU = .TRUE. LE TABLEAU U CONTIENT LA MATRICE U
!                 (MATRICE (M,N) A COLONNES ORTHOGONALES)
!                 EN SORTIE LES COLONNES SONT REORDONNEES CONFORMEMENT
!                 AUX VALEURS SINGULIERES
! IN     : MATV : LOGICAL , SCALAIRE
!                 SI MATV = .TRUE. LA MATRICE V A ETE CALCULEE
! IN/OUT : V    : REAL*8 , TABLEAU DE DIMENSION (NM,N)
!                 SI MATV = .TRUE. LE TABLEAU V CONTIENT LA MATRICE V
!                 (MATRICE CARREE D'ORDRE N ORTHOGONALE)
!                 EN SORTIE LES COLONNES SONT REORDONNEES CONFORMEMENT
!                 AUX VALEURS SINGULIERES
! IN     : EPS  : REAL*8 , SCALAIRE
!                 CRITERE DE PRECISION
! OUT    : RG   : INTEGER , SCALAIRE
!                 RANG DE LA MATRICE A
! IN/OUT : RV1  : REAL*8 , VECTEUR DE DIMENSION N
!                 VECTEUR DE TRAVAIL
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "blas/dswap.h"
    integer(kind=8) :: nm, m, n, rg
    real(kind=8) :: w(n), u(nm, n), v(nm, n), eps, rv1(n)
    aster_logical :: matu, matv
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: i, j, jmax, rgmax
    real(kind=8) :: wmax
    blas_int :: b_incx, b_incy, b_n
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   ON REORDONNE LES VALEURS SINGULIERES PAR MODULE DECROISSANT
!     SIMULTANEMENT ON EFFECTUE LES PERMUTATIONS ADEQUATES DES COLONNES
!     DES MATRICES U ET V SI ELLES ONT ETE CALCULEES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (n .gt. 1) then
        do j = 1, n-1
            jmax = j
            wmax = w(j)
            do i = j+1, n
                if (w(i) .gt. wmax) then
                    jmax = i
                    wmax = w(i)
                end if
            end do
            if (jmax .ne. j) then
                w(jmax) = w(j)
                w(j) = wmax
                if (matu) then
                    b_n = to_blas_int(m)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dswap(b_n, u(1, j), b_incx, u(1, jmax), b_incy)
                end if
                if (matv) then
                    b_n = to_blas_int(n)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dswap(b_n, v(1, j), b_incx, v(1, jmax), b_incy)
                end if
            end if
        end do
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   DETERMINATION DU RANG DE LA MATRICE A
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (w(1) .eq. 0.0d0) then
        rg = 0
    else
        rgmax = min(m, n)
        if (rgmax .gt. 1) then
            b_n = to_blas_int(rgmax)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, w(1), b_incx, rv1(1), b_incy)
            b_n = to_blas_int(rgmax)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.0d0/rv1(1), rv1(1), b_incx)
            do j = 2, rgmax
                if (rv1(j) .lt. eps) goto 40
            end do
40          continue
            rg = j-1
        else
            rg = 1
        end if
    end if
!
! --- FIN DE POSSVD.
end subroutine
