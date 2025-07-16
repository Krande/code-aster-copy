! --------------------------------------------------------------------
! Copyright (C) 2005 UCBL LYON1 - T. BARANGER     WWW.CODE-ASTER.ORG
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
subroutine hypela(fami, kpg, ksp, ndim, typmod, &
                  imate, crit, eps, option, sig, &
                  dsidep, codret)
!
    implicit none
!
#include "asterfort/hyp3ci.h"
#include "asterfort/hyp3cv.h"
#include "asterfort/hyp3di.h"
#include "asterfort/hyp3dv.h"
#include "asterfort/hypcpc.h"
#include "asterfort/hypcpd.h"
#include "asterfort/hypmat.h"
#include "asterfort/utmess.h"
#include "blas/dscal.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: imate
    real(kind=8), intent(in) :: crit(*)
    real(kind=8), intent(in) :: eps(2*ndim)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sig(6)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT HYPERELASTIQUE DE SIGNORINI
!
!     C10 (I1-3) + C01 (I2-3)+ C20 (I1-3)^2 + K/2(J-1)²
!
!     POUR LES ELEMENTS ISOPARAMETRIQUES 3D, CP, ET DP
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  TYPMOD : TYPE DE MODELISATION
! IN  CRIT   : CRITERES DE CONVERGENCE LOCAUX
!                             (1) = NB ITERATIONS MAXI A CONVERGENCE
!                                   (ITER_INTE_MAXI == ITECREL)
!                             (2) = TYPE DE JACOBIEN A T+DT
!                                   (TYPE_MATR_COMP == MACOMP)
!                                   0 = EN VITESSE     >SYMETRIQUE
!                                   1 = EN INCREMENTAL >NON-SYMETRIQUE
!                             (3) = VALEUR TOLERANCE DE CONVERGENCE
!                                    (RESI_INTE == RESCREL)
! IN  IMATE  : ADRESSE DU MATERIAU CODE
! IN  FAMI   : FAMILLE DE POINTS DE GAUSS
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  KSP    : NUMERO DU SOUS-POINT DE GAUSS
! IN  EPS    : DEFORMATION (SI C_PLAN EPS(3) EST EN FAIT CALCULE)
! OUT SIG    : CONTRAINTES
! OUT DSIDEP : MATRICE TANGENTE
! OUT CODRET : CODE RETOUR ERREUR INTEGRATION (1 SI PROBLEME, 0 SINON)
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, l, m
    real(kind=8) :: c11, c22, c12, c33, c13, c23
    real(kind=8) :: epstot(6)
    real(kind=8) :: cvol(6, 6), ciso(6, 6)
    real(kind=8) :: svol(6), siso(6)
    real(kind=8) :: c10, c01, c20, k
    integer(kind=8) :: nitmax
    real(kind=8) :: epsi
    character(len=1) :: poum
    blas_int :: b_incx, b_n
!
! ----------------------------------------------------------------------
!
    dsidep(:, :) = 0.d0
    poum = merge('-', '+', option(1:9) .eq. 'RIGI_MECA')
!
! --- LECTURE DES CARACTERISTIQUES MATERIAU
!
    call hypmat(fami, kpg, ksp, poum, imate, &
                c10, c01, c20, k)
!
! --- A PRIORI ON A CONVERGE
!
    codret = 0
!
! --- PRE-TRAITEMENT DES DEFORMATIONS (PAS DE NOTATION DE VOIGT)
!
    epstot = 0
    do i = 1, 3
        epstot(i) = eps(i)
    end do
    do i = 4, 2*ndim
        epstot(i) = eps(i)/sqrt(2.d0)
    end do
!
! --- CALCUL CONTRAINTES ET MATRICE TANGENTE
!
    if (typmod(1) .eq. '3D' .or. typmod(1) .eq. '3D_SI' .or. typmod(1) (1:6) .eq. 'D_PLAN') then
! --- CALCUL DES ELONGATIONS
        c11 = 2.d0*epstot(1)+1.d0
        c12 = 2.d0*epstot(4)
        c22 = 2.d0*epstot(2)+1.d0
        c33 = 2.d0*epstot(3)+1.d0
        c13 = 2.d0*epstot(5)
        c23 = 2.d0*epstot(6)
! --- CALCUL DES CONTRAINTES ISOTROPIQUES
        call hyp3ci(c11, c22, c33, c12, c13, &
                    c23, c10, c01, c20, siso, &
                    codret)
        if (codret .eq. 1) then
            goto 99
        end if
! --- CALCUL DES CONTRAINTES VOLUMIQUES
        call hyp3cv(c11, c22, c33, c12, c13, &
                    c23, k, svol, codret)
        if (codret .eq. 1) then
            goto 99
        end if
! --- CALCUL DE LA MATRICE TANGENTE (PARTIE ISOTROPIQUE)
        call hyp3di(c11, c22, c33, c12, c13, &
                    c23, c10, c01, c20, ciso, &
                    codret)
        if (codret .eq. 1) then
            goto 99
        end if
! --- CALCUL DE LA MATRICE TANGENTE (PARTIE VOLUMIQUE)
        call hyp3dv(c11, c22, c33, c12, c13, &
                    c23, k, cvol, codret)
        if (codret .eq. 1) then
            goto 99
        end if
!
! --- ASSEMBLAGE VOLUMIQUE/ISOTROPIQUE
! --- ON CORRIGE A CE NIVEAU LES TERMES LIES AU CISAILLEMENT
! --- A TERME IL FAUDRA RE-ECRIRE LES ROUTINES D'INTEGRATION
        do i = 1, 3
            sig(i) = siso(i)+svol(i)
            sig(3+i) = (siso(3+i)+svol(3+i))/2.d0
            do j = 1, 3
                dsidep(i, j) = ciso(i, j)+cvol(i, j)
                dsidep(3+i, j) = (ciso(3+i, j)+cvol(3+i, j))/sqrt(2.d0)
                dsidep(i, 3+j) = (ciso(i, 3+j)+cvol(i, 3+j))/sqrt(2.d0)
                dsidep(3+i, 3+j) = (ciso(3+i, 3+j)+cvol(3+i, 3+j))/2.d0
            end do
        end do
    else if (typmod(1) (1:6) .eq. 'C_PLAN') then
        epsi = abs(crit(3))
        nitmax = abs(nint(crit(1)))
!
! --- CALCUL DES ELONGATIONS
        c11 = 2.d0*epstot(1)+1.d0
        c12 = 2.d0*epstot(4)
        c22 = 2.d0*epstot(2)+1.d0
        c33 = 1.d0
! --- CALCUL DES CONTRAINTES
        call hypcpc(c11, c22, c33, c12, k, &
                    c10, c01, c20, nitmax, epsi, &
                    sig, codret)
        if (codret .eq. 1) then
            goto 99
        end if
! --- CALCUL DE LA MATRICE TANGENTE
        call hypcpd(c11, c22, c33, c12, k, &
                    c10, c01, c20, dsidep, codret)
        if (codret .eq. 1) then
            goto 99
        end if
! --- ON CORRIGE A CE NIVEAU LES TERMES LIES AU CISAILLEMENT
! --- A TERME IL FAUDRA RE-ECRIRE LES ROUTINES D'INTEGRATION
        do i = 1, 3
            sig(3+i) = sig(3+i)/2.d0
            do j = 1, 3
                dsidep(3+i, j) = dsidep(3+i, j)/sqrt(2.d0)
                dsidep(i, 3+j) = dsidep(i, 3+j)/sqrt(2.d0)
                dsidep(3+i, 3+j) = dsidep(3+i, 3+j)/2.d0
            end do
        end do
! --- PRISE EN COMPTE DE L'HYPOTHESE DE CONTRAINTES PLANES DANS DSIDEP
        do m = 1, 2*ndim
            if (m .ne. 3) then
                do l = 1, 2*ndim
                    if (l .ne. 3) then
                        dsidep(m, l) = dsidep(m, l)-1.d0/dsidep(3, 3)*dsidep(m, 3)*dsidep(3, l)
                    end if
                end do
            end if
        end do
    else
        call utmess('F', 'ELASHYPER_97', sk=typmod(1))
    end if
!
! --- POST-TRAITEMENT DES CONTRAINTES (PAS DE NOTATION DE VOIGT)
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, sqrt(2.d0), sig(4), b_incx)
!
99  continue
end subroutine
