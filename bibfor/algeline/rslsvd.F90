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
subroutine rslsvd(nm, m, n, a, w, &
                  u, v, nb, b, eps, &
                  ierr, rvnm)
    implicit none
!
! DESCRIPTION :   RESOLUTION D'UN SYSTEME LINEAIRE RECTANGULAIRE (M,N)
! -----------                          A(X) = B
!                 AU SENS DES MOINDRES CARRES, PAR DECOMPOSITION AUX
!                 VALEURS SINGULIERES DE LA MATRICE A
!
! IN     : NM   : INTEGER , SCALAIRE
!                 PREMIERE DIMENSION DES TABLEAUX A, U ET V, DECLAREE
!                 DANS L'APPELANT, NM >= MAX(M,N)
! IN     : M    : INTEGER , SCALAIRE
!                 NOMBRE DE LIGNES DES MATRICES A ET U
!                 NOMBRE DE COLONNES DE U
! IN     : N    : INTEGER , SCALAIRE
!                 NOMBRE DE COLONNES DE A
!                  = ORDRE DE LA MATRICE V
! IN     : A    : REAL*8 , TABLEAU DE DIMENSION(NM,N)
!                 CONTIENT LA MATRICE RECTANGULAIRE A
!                 LE CONTENU EST INCHANGE EN SORTIE
! OUT    : W    : REAL*8 , VECTEUR DE DIMENSION N
!                 CONTIENT LES N VALEURS SINGULIERES DE A, ORDONNEES
!                 PAR MODULE DECROISSANT
! OUT    : U    : REAL*8 , TABLEAU DE DIMENSION (NM,M)
!                 CONTIENT LA MATRICE U, MATRICE (M,M) A COLONNES
!                 ORTHOGONALES
! OUT    : V    : REAL*8 , TABLEAU DE DIMENSION (NM,N)
!                 CONTIENT LA MATRICE V, MATRICE CARREE D'ORDRE N
!                 ORTHOGONALE
! IN     : NB   : INTEGER , SCALAIRE
!                 NOMBRE DE SECONDS MEMBRES, I.E. NOMBRE DE COLONNES DE
!                 LA MATRICE B
! IN/OUT : B    : REAL*8 , TABLEAU DE DIMENSION (NM,NB)
!                 EN ENTREE : LES VECTEURS COLONNES CORRESPONDENT AUX
!                 SECONDS MEMBRES (VECTEURS D'ORDRE M)
!                 EN SORTIE : LES VECTEURS COLONNES CORRESPONDENT AUX
!                 SOLUTIONS (VECTEURS D'ORDRE N)
! IN     : EPS  : REAL*8 , SCALAIRE
!                 CRITERE DE PRECISION
! OUT    : IERR : INTEGER , SCALAIRE , CODE RETOUR
!                 IERR = 0  : OK
!                 IERR /= 0 : PB LORS DE LA DECOMPOSITION
! IN/OUT : RVNM : REAL*8 , VECTEUR DE DIMENSION NM
!                 VECTEUR DE TRAVAIL
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "asterfort/calsvd.h"
#include "asterfort/possvd.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
    integer(kind=8) :: nm, m, n, nb, ierr
    real(kind=8) :: a(nm, n), w(n), u(nm, m), v(nm, n), b(nm, nb), eps, rvnm(nm)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: ib, j, rg
    real(kind=8) :: alphaj
    aster_logical :: matuv
    blas_int :: b_incx, b_incy, b_n
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   DECOMPOSITION AUX VALEURS SINGULIERES DE LA MATRICE A
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    matuv = .true.
    call calsvd(nm, m, n, a, w, &
                matuv, u, matuv, v, ierr)
    if (ierr .ne. 0) goto 999
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   POST-TRAITEMENTS AU CALCUL DE LA DECOMPOSITION AUX VALEURS
!     SINGULIERES DE LA MATRICE A :
!     LES VALEURS SINGULIERES SONT REORDONNEES PAR MODULE DECROISSANT
!     SIMULTANEMENT ON EFFECTUE LES PERMUTATIONS ADEQUATES DES COLONNES
!     DES MATRICES U ET V
!     DETERMINATION DU RANG DE LA MATRICE A
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    call possvd(nm, m, n, w, matuv, &
                u, matuv, v, eps, rg, &
                rvnm)
    if (rg .eq. 0) then
        ierr = -1
        goto 999
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   CALCUL DES SOLUTIONS DU SYSTEME
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    do ib = 1, nb
        call r8inir(n, 0.0d0, rvnm(1), 1)
        do j = 1, rg
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            alphaj = ddot(b_n, u(1, j), b_incx, b(1, ib), b_incy)/w(j)
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, alphaj, v(1, j), b_incx, rvnm(1), &
                       b_incy)
        end do
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, rvnm(1), b_incx, b(1, ib), b_incy)
    end do
!
999 continue
!
! --- FIN DE RSLSVD.
end subroutine
