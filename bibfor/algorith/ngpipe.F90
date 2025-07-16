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
subroutine ngpipe(typilo, npg, neps, nddl, b, &
                  ni2ldc, typmod, mat, compor, lgpg, &
                  ddlm, sigm, vim, ddld, ddl0, &
                  ddl1, tau, etamin, etamax, copilo)
!
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/pil000.h"
#include "blas/dgemv.h"
    character(len=8) :: typmod(*)
    character(len=16) :: typilo, compor(*)
!
    integer(kind=8) :: npg, neps, nddl, mat, lgpg
    real(kind=8) :: ddlm(nddl), ddld(nddl), ddl0(nddl), ddl1(nddl)
    real(kind=8) :: sigm(neps, npg), vim(lgpg, npg), tau
    real(kind=8) :: copilo(5, npg), etamin, etamax
    real(kind=8) :: b(neps, npg, nddl), ni2ldc(neps, npg)
!.......................................................................
!
!     BUT:  CALCUL  DES COEFFICIENTS DE PILOTAGE POUR PRED_ELAS
!.......................................................................
! IN  TYPILO : MODE DE PILOTAGE: 'DEFORMATION', 'PRED_ELAS'
! IN  NPG    : NOMBRE DE POINTS DE GAUSS
! IN  NEPS   : NOMBRE DE COMPOSANTES DE DEFORMATIONS / CONTRAINTES
! IN  NDDL   : NOMBRE DE DDL DANS L'ELEMENT
! IN  B      : MATRICE CINEMATIQUE
! IN  ni2ldc : CONVERSION CONTRAINTE --> AVEC RACINE DE DEUX
! IN  TYPMOD : TYPE DE MODELISATION
! IN  MAT    : MATERIAU CODE
! IN  COMPOR : COMPORTEMENT
! IN  LGPG   : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!             CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  DDLM   : DDL U,ALPHA,MU EN T-
! IN  SIGM   : CONTRAINTES DE CAUCHY EN T- (INUTILE)
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  DDLD   : INCREMENT DE DDL U,ALPHA,MU A L'ITERATION NEWTON COURANTE
! IN  DDL0   : CORRECTION DE DDL U,ALPHA,MU POUR FORCES FIXES
! IN  DDL1   : CORRECTION DE DDL U,ALPHA,MU POUR FORCES PILOTEES
! OUT COPILO : COEFFICIENTS A0 ET A1 POUR CHAQUE POINT DE GAUSS
! ----------------------------------------------------------------------
    integer(kind=8) :: g, nepg
    real(kind=8) :: sigmam(neps, npg)
    real(kind=8) :: epsm(neps, npg), epsd_pilo(neps, npg)
    real(kind=8) :: epsd_cste(neps, npg)
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
! ----------------------------------------------------------------------
!
!
! -- INITIALISATION
!
    ASSERT(compor(3) .eq. 'PETIT')
    copilo = r8vide()
    nepg = neps*npg
!
!
! -- DEFORMATIONS
!
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddlm, b_incx, 0.d0, epsm, &
               b_incy)
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddld, b_incx, 0.d0, epsd_cste, &
               b_incy)
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddl0, b_incx, 1.d0, epsd_cste, &
               b_incy)
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddl1, b_incx, 0.d0, epsd_pilo, &
               b_incy)
!
!
! -- PRETRAITEMENT SI NECESSAIRE
    if (typilo .eq. 'PRED_ELAS') sigmam = sigm*ni2ldc

!
! -- TRAITEMENT DE CHAQUE POINT DE GAUSS
!
    do g = 1, npg
        call pil000(typilo, compor, neps, tau, mat, &
                    vim(:, g), sigmam(1, g), epsm(1, g), epsd_cste(1, g), epsd_pilo(1, g), &
                    typmod, etamin, etamax, copilo(1, g))
    end do
!
end subroutine
