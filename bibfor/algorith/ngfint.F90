! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504,W1306
!
subroutine ngfint(option, typmod, ndim, nddl, neps, &
                  npg, w, b, compor, fami, &
                  mat, angmas, lgpg, carcri, instam, &
                  instap, ddlm, ddld, ni2ldc, sigmam, &
                  vim, sigmap, vip, fint, matr, &
                  lMatr, lVect, lSigm, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codere.h"
#include "asterfort/nmcomp.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
!
    character(len=8) :: typmod(2)
    character(len=*) :: fami
    character(len=16) :: option, compor(COMPOR_SIZE)
    integer :: ndim, nddl, neps, npg, mat, lgpg
    real(kind=8) :: w(neps, npg), ni2ldc(neps, npg), b(neps, npg, nddl)
    real(kind=8) :: angmas(3), carcri(*), instam, instap
    real(kind=8) :: ddlm(nddl), ddld(nddl)
    real(kind=8) :: sigmam(neps, npg), sigmap(neps, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), matr(nddl, nddl), fint(nddl)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
    integer, intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     RAPH_MECA, RIGI_MECA_* ET FULL_MECA_*
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODEELISATION                              (LDC)
! IN  NDIM    : DIMENSION DE L'ESPACE                              (LDC)
! IN  NDDL    : NOMBRE DE DEGRES DE LIBERTE
! IN  NEPS    : NOMBRE DE COMPOSANTES DE DEFORMATION ET CONTRAINTE
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  W       : POIDS DES POINTS DE GAUSS
! IN  B       : MATRICE CINEMATIQUE : DEFORMATION = B.DDL
! IN  COMPOR  : COMPORTEMENT                                       (LDC)
! IN  MAT     : MATERIAU CODE                                      (LDC)
! IN  ANGMAS  : ANGLE DU REPERE LOCAL                              (LDC)
! IN  LGPG    : LONGUEUR DU TABLEAU DES VARIABLES INTERNES
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX                     (LDC)
! IN  INSTAM  : INSTANT PRECEDENT                                  (LDC)
! IN  INSTAP  : INSTANT DE CALCUL                                  (LDC)
! IN  DDLM    : DDL A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DDL
! IN  LI2LDC  : CONVERSION CONTRAINTE STOCKEE -> CONTRAINTE LDC (RAC2)
! IN  SIGMAM  : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGMAP  : CONTRAINTES DE CAUCHY (RAPH_MECA   ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA   ET FULL_MECA_*)
! OUT FINT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! OUT MATR    : MATRICE DE RIGIDITE   (RIGI_MECA_* ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: ksp = 1
    integer :: nepg, g, i, cod(npg)
    real(kind=8) :: sigm(neps, npg), sigp(neps, npg)
    real(kind=8) :: epsm(neps, npg), epsd(neps, npg)
    real(kind=8) :: dsidep(neps, neps, npg)
    real(kind=8) :: ktgb(0:neps*npg*nddl-1)
    type(Behaviour_Integ) :: BEHinteg
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    blas_int :: b_incx, b_incy
!
! --------------------------------------------------------------------------------------------------
!
    nepg = neps*npg
    codret = 0
    if (lMatr) dsidep = 0.d0
    if (lSigm) sigp = 0.d0
    cod = 0

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, carcri, &
                              instam, instap, &
                              fami, mat, &
                              BEHinteg)
!
! - CALCUL DES DEFORMATIONS GENERALISEES
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
               b_lda, ddld, b_incx, 0.d0, epsd, &
               b_incy)
!
! - CALCUL DE LA LOI DE COMPORTEMENT
!
!    FORMAT LDC DES CONTRAINTES (AVEC RAC2)
    sigm = sigmam*ni2ldc
!
! - LOI DE COMPORTEMENT EN CHAQUE POINT DE GAUSS
    do g = 1, npg

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(g, ksp, BEHinteg)

! ----- Integrator
        call nmcomp(BEHinteg, &
                    fami, g, ksp, ndim, typmod, &
                    mat, compor, carcri, instam, instap, &
                    neps, epsm(:, g), epsd(:, g), neps, sigm(:, g), &
                    vim(1, g), option, angmas, &
                    sigp(:, g), vip(1, g), neps*neps, dsidep(:, :, g), cod(g))
        if (cod(g) .eq. 1) goto 900
    end do
!
!    FORMAT RESULTAT DES CONTRAINTES (SANS RAC2)
    if (lSigm) sigmap = sigp/ni2ldc
!
! - FORCE INTERIEURE
!
    if (lVect) then
!      PRISE EN CHARGE DU POIDS DU POINT DE GAUSS
        sigp = sigp*w
!      FINT = SOMME(G) WG.BT.SIGMA
        b_lda = to_blas_int(nepg)
        b_m = to_blas_int(nepg)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, 1.d0, b, &
                   b_lda, sigp, b_incx, 0.d0, fint, &
                   b_incy)
    end if
!
! - CALCUL DE LA MATRICE DE RIGIDITE (STOCKAGE PAR LIGNES SUCCESSIVES)
!
    if (lMatr) then
!      PRISE EN CHARGE DU POIDS DU POINT DE GAUSS  WG.DSIDEP
        do i = 1, neps
            dsidep(:, i, :) = dsidep(:, i, :)*w
        end do
!      CALCUL DES PRODUITS INTERMEDIAIRES (WG.DSIDEP).B POUR CHAQUE G
        do g = 1, npg
            b_ldc = to_blas_int(nepg)
            b_ldb = to_blas_int(nepg)
            b_lda = to_blas_int(neps)
            b_m = to_blas_int(neps)
            b_n = to_blas_int(nddl)
            b_k = to_blas_int(neps)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, dsidep(1, 1, g), b_lda, b(1, g, 1), b_ldb, &
                       0.d0, ktgb((g-1)*neps), b_ldc)
        end do
!      CALCUL DU PRODUIT FINAL SOMME(G) BT. ((WG.DSIDEP).B)  TRANSPOSE
        b_ldc = to_blas_int(nddl)
        b_ldb = to_blas_int(nepg)
        b_lda = to_blas_int(nepg)
        b_m = to_blas_int(nddl)
        b_n = to_blas_int(nddl)
        b_k = to_blas_int(nepg)
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   1.d0, ktgb, b_lda, b, b_ldb, &
                   0.d0, matr, b_ldc)
    end if
!
! - SYNTHESE DU CODE RETOUR
900 continue
    if (lSigm) call codere(cod, npg, codret)
!
end subroutine
