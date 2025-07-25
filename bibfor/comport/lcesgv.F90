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
subroutine lcesgv(fami, kpg, ksp, ndim, neps, &
                  typmod, option, mat, lccrma, lcesga, &
                  epsm, deps, vim, itemax, precvg, &
                  sig, vip, dsidep, iret)
    implicit none
#include "asterf_types.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "asterfort/assert.h"
#include "asterfort/lcesel.h"
#include "asterfort/lcesma.h"
#include "asterfort/lcesrf.h"
#include "asterfort/lcesvf.h"
#include "asterfort/lcgrad.h"
#include "asterc/r8prem.h"
    interface
        subroutine lccrma(mat, fami, kpg, ksp, poum)
            integer(kind=8), intent(in) :: mat, kpg, ksp
            character(len=1), intent(in) :: poum
            character(len=*), intent(in) :: fami
        end subroutine lccrma
!
        subroutine lcesga(mode, eps, gameps, dgamde, itemax, &
                          precvg, iret)
            integer(kind=8), intent(in) :: mode, itemax
            real(kind=8), intent(in) :: eps(6), precvg
            integer(kind=8), intent(out) :: iret
            real(kind=8), intent(out) :: gameps, dgamde(6)
        end subroutine lcesga
    end interface
!
    character(len=8) :: typmod(*)
    character(len=16) :: option
    character(len=*) :: fami
    integer(kind=8) :: ndim, neps, mat, iret, kpg, ksp, itemax
    real(kind=8) :: epsm(neps), deps(neps), vim(*), precvg
    real(kind=8) :: vip(*), sig(neps), dsidep(neps, neps)
! --------------------------------------------------------------------------------------------------
!           ENDOMMAGEMENT FRAGILE A GRADIENT DE VARIABLE INTERNE :
!                       ENDO_SCALAIRE AVEC GRAD_VARI
! --------------------------------------------------------------------------------------------------
! IN  NEPS    DIMENSION DES DEFORMATIONS GENERALISEES
! IN  TYPMOD  TYPE DE MODELISATION
! IN  OPTION  OPTION DE CALCUL
!               RIGI_MECA_TANG, RIGI_MECA_ELAS
!               RAPH_MECA
!               FULL_MECA, FULL_MECA_ELAS
! IN  MAT     NATURE DU MATERIAU
! IN  LCCRMA  ROUTINE POUR LECTURE MATERIAU SPECIFIQUE AU CRITERE
! IN  LCGVGA  ROUTINE POUR CALCUL GAMMA(EPS) SPECIFIQUE AU CRITERE
! IN  EPSM    CHAMP DE DEFORMATION EN T- ET PHIM=EPSM(7)
! IN  DEPS    INCREMENT DU CHAMP DE DEFORMATION ET DPHI=DEPS(7)
! IN  VIM     VARIABLES INTERNES EN T-
! IN  NONLOC  INUTILISE
! IN  ITEMAX  NBR MAXI D'ITERATIONS POUR RESOLUTION EQUATION SCALAIRE
! IN  PRECVG  CRITERE DE CVG : 10 A 100 FOIS PLUS FIN QUE RESI_REFE_RELA
! OUT VIP     DENSITE DE FISSURATION
! OUT SIG     CONTRAINTE
! OUT DSIDEP  MATRICE TANGENTE
! OUT IRET    CODE RETOUR (0=OK, 1=ECHEC CVG)
! --------------------------------------------------------------------------------------------------
    real(kind=8), dimension(6), parameter :: kr = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
! --------------------------------------------------------------------------------------------------
    aster_logical :: cplan, rigi, resi, elas
    integer(kind=8) :: ndimsi, ij, kl, etat
    real(kind=8) :: phi, lag, apg, grad(ndim)
    real(kind=8) :: coplan, cor33, vplan(6), eps(6), sigel(6), treps
    real(kind=8) :: sigma(6), a, drda, drdae, drdas, gel, gsat, ktg(6, 6, 4)
    real(kind=8) :: ra, fd, d2rda2, dgda, gameps, dgamde(6), coefg
    real(kind=8) :: yng, nrmela, sigref, a0, d1a0, prece, preca, precga
    real(kind=8) :: sigela(6), sigelu(6), dsade(6, 6), dsude(6, 6)
    character(len=1) :: poum
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
    common/lcee/lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: pk, pm, pp, pq
    blas_int :: b_incx, b_incy, b_n
    common/lces/pk, pm, pp, pq
! --------------------------------------------------------------------------------------------------
!
!
! --------------------------------------------------------------------------------------------------
!                          INITIALISATIONS
! --------------------------------------------------------------------------------------------------
!
!
! -- OPTIONS DE CALCUL
!
    cplan = typmod(1) .eq. 'C_PLAN  '
    elas = option(11:14) .eq. 'ELAS'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    resi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH'
    eps(5) = 0
    eps(6) = 0
    ndimsi = 2*ndim
    iret = 0
    poum = merge('+', '-', resi)
!
!
! -- LECTURE DES CARACTERISTIQUES MATERIAU
!
    call lcesma(mat, fami, kpg, ksp, poum, &
                lccrma)
    ASSERT(gamma .eq. 0.d0 .or. .not. cplan)
!
!
! -- DEFORMATIONS COURANTES
!
!    DEFORMATION, ENDOMMAGEMENT, LAGRANGE ET GRADIENT
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, eps, b_incy)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm(ndimsi+3), b_incx, grad, b_incy)
    apg = epsm(ndimsi+1)
    lag = epsm(ndimsi+2)
!
    if (resi) then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, eps, &
                   b_incy)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps(ndimsi+3), b_incx, grad, &
                   b_incy)
        apg = apg+deps(ndimsi+1)
        lag = lag+deps(ndimsi+2)
    end if
!
    phi = lag+pr*apg
!
!    DEFORMATIONS MECANIQUES
    eps(1) = eps(1)-epsth
    eps(2) = eps(2)-epsth
    eps(3) = eps(3)-epsth
!
!    CONTRAINTES PLANES
    if (cplan) then
        coplan = -lambda/(lambda+deuxmu)
        eps(3) = coplan*(eps(1)+eps(2))
    end if
!
!
!  CONTRAINTE ELASTIQUE
    treps = eps(1)+eps(2)+eps(3)
    sigel = lambda*treps*kr+deuxmu*eps
!
!
!  SEUIL DES CRITERES DE CONVERGENCE
    yng = deuxmu*(3*lambda+deuxmu)/(2*lambda+deuxmu)
    sigref = sqrt(2*yng*pk/pm)
    nrmela = sqrt(dot_product(sigel, sigel))
    a0 = vim(1)
    d1a0 = abs(lcesvf(1, a0))
    prece = precvg*sigref/yng
    preca = precvg*minval((/1.d0, pk/pr, sigref/max(r8prem(), d1a0*nrmela)/))
    precga = pr*preca/max(r8prem(), d1a0)
!
!
!
! -- PSEUDO-ENERGIE DE DEFORMATION ET CONTRAINTE ELASTIQUE
!
!
    call lcesga(0, eps, gameps, dgamde, itemax, &
                precga, iret)
    if (iret .ne. 0) goto 999
    if (resi) then
        call lcesel(eps, rigi, elas, prece, sigela, &
                    sigelu, dsade, dsude)
    else
        call lcesel(vim(4:9), rigi, elas, prece, sigela, &
                    sigelu, dsade, dsude)
    end if
!
!
!
! --------------------------------------------------------------------------------------------------
!                     CALCUL DE L'ENDOMMAGEMENT
! --------------------------------------------------------------------------------------------------
!
    a = vim(1)
    etat = nint(vim(2))
!
    if (.not. resi) goto 500
!
!
!    PAS DE CALCUL D'ENDOMMAGEMENT POUR UN POINT SATURE
    if (etat .eq. 2) goto 200
!
!
!    ESTIMATION DU CRITERE
!
!    PREDICTION ELASTIQUE
!
    drdae = lcesvf(1, a)
    gel = drdae*gameps+pk-phi+pr*a
!
    if (gel .ge. 0) then
        etat = 0
        goto 200
    end if
!
!
!    PREDICTION SATUREE
!
    drdas = lcesvf(1, 1.d0)
    gsat = drdas*gameps+pk-phi+pr
    if (gsat .le. 0) then
        etat = 2
        a = 1.d0
        goto 200
    end if
!
!
!    RESOLUTION DE L'EQUATION G(A)=0
    etat = 1
    a = lcesrf(a, gameps, pr, pk-phi, preca, itemax, iret)
    if (iret .ne. 0) goto 999
!
!    PROJECTION DE A+ ENTRE A- ET 1.D0
    if (a .le. vim(1)) then
        a = vim(1)
    else if (a .gt. 1.d0) then
        etat = 2
        a = 1.d0
    end if
!
!
!    STOCKAGE DES CONTRAINTES ET DES VARIABLES INTERNES
!
200 continue
!
    ra = lcesvf(0, a)
    sigma = ra*sigela+sigelu
!
    vip(1) = a
    vip(2) = etat
    vip(3) = 1.d0-ra
    vip(4:9) = eps
!
500 continue
!
!
! --------------------------------------------------------------------------------------------------
!                     CALCUL DES MATRICES TANGENTES
! --------------------------------------------------------------------------------------------------
!
    if (.not. rigi) goto 800
!
    ktg = 0
!
!
! -- CONTRIBUTION ELASTIQUE
!
    ra = lcesvf(0, a)
    fd = max(ra, rigmin)
    ktg(:, :, 1) = fd*dsade+dsude
!
!
!
! -- CORRECTION DISSIPATIVE
    if (etat .eq. 1 .and. .not. elas) then
!
        call lcesga(1, eps, gameps, dgamde, itemax, &
                    precga, iret)
        drda = lcesvf(1, a)
        d2rda2 = lcesvf(2, a)
        dgda = d2rda2*gameps+pr
        coefg = drda**2/dgda
!
        do ij = 1, ndimsi
            do kl = 1, ndimsi
                ktg(ij, kl, 1) = ktg(ij, kl, 1)-coefg*sigela(ij)*dgamde(kl)
            end do
            ktg(ij, 1, 2) = drda/dgda*sigela(ij)
            ktg(ij, 1, 3) = -drda/dgda*dgamde(ij)
        end do
        ktg(1, 1, 4) = 1/dgda
!
    end if
!
!
! -- CORRECTION POUR LES CONTRAINTES PLANES
!
    if (cplan) then
!
        cor33 = coplan**2*ktg(3, 3, 1)
        do ij = 1, ndimsi
            vplan(ij) = coplan*ktg(ij, 3, 1)
        end do
        do ij = 1, ndimsi
            ktg(ij, 1, 1) = ktg(ij, 1, 1)+vplan(ij)
            ktg(ij, 2, 1) = ktg(ij, 2, 1)+vplan(ij)
            ktg(1, ij, 1) = ktg(1, ij, 1)+vplan(ij)
            ktg(2, ij, 1) = ktg(2, ij, 1)+vplan(ij)
        end do
        ktg(1, 1, 1) = ktg(1, 1, 1)+cor33
        ktg(1, 2, 1) = ktg(1, 2, 1)+cor33
        ktg(2, 1, 1) = ktg(2, 1, 1)+cor33
        ktg(2, 2, 1) = ktg(2, 2, 1)+cor33
!
        ktg(1, 1, 2) = ktg(1, 1, 2)+coplan*ktg(3, 1, 2)
        ktg(2, 1, 2) = ktg(2, 1, 2)+coplan*ktg(3, 1, 2)
        ktg(1, 1, 3) = ktg(1, 1, 3)+coplan*ktg(3, 1, 3)
        ktg(2, 1, 3) = ktg(2, 1, 3)+coplan*ktg(3, 1, 3)
!
    end if
!
!
! -- PRISE EN CHARGE DES TERMES DU LAGRANGIEN AUGMENTE
!
800 continue
    call lcgrad(resi, rigi, sigma(1:ndimsi), apg, lag, &
                grad, a, pr, pc, ktg(1:ndimsi, 1:ndimsi, 1), &
                ktg(1:ndimsi, 1, 2), ktg(1:ndimsi, 1, 3), ktg(1, 1, 4), sig, dsidep)
!
!
999 continue
end subroutine
