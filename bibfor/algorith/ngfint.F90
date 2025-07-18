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
! aslint: disable=W1504,W1306
!
subroutine ngfint(option, typmod, ndim, nddl, neps, &
                  npg, w, b, compor, fami, &
                  mat, angmas, lgpg, carcri, instam, &
                  instap, ddlm, ddld, ni2ldc, sigmam, &
                  vim, sigmap, vip, fint, matsym, matuu, matns, &
                  lMatr, lVect, lSigm, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/codere.h"
#include "asterfort/nmcomp.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
!
    character(len=8) :: typmod(2)
    character(len=*) :: fami
    character(len=16) :: option, compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, nddl, neps, npg, mat, lgpg
    real(kind=8) :: w(neps, npg), ni2ldc(neps, npg), b(neps, npg, nddl)
    real(kind=8) :: angmas(3), carcri(*), instam, instap
    real(kind=8) :: ddlm(nddl), ddld(nddl)
    real(kind=8) :: sigmam(neps, npg), sigmap(neps, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), fint(nddl)
    real(kind=8), intent(out):: matuu((nddl*(nddl+1))/2)
    real(kind=8), intent(out), target :: matns(nddl, nddl)
    aster_logical, intent(in) :: matsym, lMatr, lVect, lSigm
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     RAPH_MECA, RIGI_MECA_* ET FULL_MECA_*
!
! --------------------------------------------------------------------------------------------------
! in  option  : option de calcul
! in  typmod  : type de modeelisation                              (ldc)
! in  ndim    : dimension de l'espace                              (ldc)
! in  nddl    : nombre de degres de liberte
! in  neps    : nombre de composantes de deformation et contrainte
! in  npg     : nombre de points de gauss
! in  w       : poids des points de gauss
! in  b       : matrice cinematique : deformation = b.ddl
! in  compor  : comportement                                       (ldc)
! in  mat     : materiau code                                      (ldc)
! in  angmas  : angle du repere local                              (ldc)
! in  lgpg    : longueur du tableau des variables internes
! in  crit    : criteres de convergence locaux                     (ldc)
! in  instam  : instant precedent                                  (ldc)
! in  instap  : instant de calcul                                  (ldc)
! in  ddlm    : ddl a l'instant precedent
! in  ddld    : increment des ddl
! in  li2ldc  : conversion contrainte stockee -> contrainte ldc (rac2)
! in  sigmam  : contraintes a l'instant precedent
! in  vim     : variables internes a l'instant precedent
! out sigmap  : contraintes de cauchy (raph_meca   et full_meca_*)
! out vip     : variables internes    (raph_meca   et full_meca_*)
! out fint    : forces interieures    (raph_meca   et full_meca_*)
! in  matsym  : doit-on remplir la matrice symétrique ou la matrice non-symétrique
! out matuu   : matrice de rigidite symetrique   (rigi_meca_* et full_meca_*)
! out matns   : matrice de rigidite non symetrique (rigi_meca_* et full_meca_*)
! out codret  : code retour
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: nepg, g, i, j, cod(npg)
    real(kind=8) :: sigm(neps, npg), sigp(neps, npg)
    real(kind=8) :: epsm(neps, npg), epsd(neps, npg)
    real(kind=8) :: dsidep(neps, neps, npg)
    real(kind=8) :: ktgb(0:neps*npg*nddl-1)
    real(kind=8), pointer, dimension(:, :) :: ktan_t => null()
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

        ! Matrice symétrique ou non
        if (matsym) then
            allocate (ktan_t(nddl, nddl))
        else
            ktan_t => matns
        end if

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
                   0.d0, ktan_t, b_ldc)

        ! Stockage de la matrice symetrique
        if (matsym) then
            forall (i=1:nddl, j=1:nddl, i .ge. j) matuu((i*(i-1))/2+j) = ktan_t(j, i)
            deallocate (ktan_t)
        end if
    end if
!
! - SYNTHESE DU CODE RETOUR
900 continue
    if (lSigm) call codere(cod, npg, codret)
!
end subroutine
