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
subroutine gdclci(fm, df, em)
!
!
    implicit none
!
#include "blas/ddot.h"
    real(kind=8) :: fm(3, 3), df(3, 3), em(6)
! ----------------------------------------------------------------------
!       INTEGRATION DES LOIS EN GRANDES DEFORMATIONS CANO-LORENTZ
!                  CALCUL DES ELEMENTS CINEMATIQUES
! ----------------------------------------------------------------------
! IN  FM    DEFORMATION AU DEBUT DU PAS DE TEMPS
! IN  DF    INCREMENT DE DEFORMATION PENDANT LE PAS DE TEMPS
! IN  EM    DEFORMATION ELASTIQUE (-) AU DEBUT DU PAS DE TEMPS
! ----------------------------------------------------------------------
!  COMMON GRANDES DEFORMATIONS CANO-LORENTZ
!
    integer(kind=8) :: ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6)
    real(kind=8) :: lambda, mu, deuxmu, unk, troisk, cother
    real(kind=8) :: jm, dj, jp, djdf(3, 3)
    real(kind=8) :: etr(6), dvetr(6), eqetr, tretr, detrdf(6, 3, 3)
    real(kind=8) :: dtaude(6, 6)
!
    common/gdclc/&
     &          ind1, ind2, kr, rac2, rc,&
     &          lambda, mu, deuxmu, unk, troisk, cother,&
     &          jm, dj, jp, djdf,&
     &          etr, dvetr, eqetr, tretr, detrdf,&
     &          dtaude
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: ij, kl, i, j, k, l
    real(kind=8) :: bem(6), pdf(6, 6), betr(6)
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
!
!
!  CALCUL DES JACOBIENS
! ----------------------
!
    jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2))-fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(&
         &3, 2))+fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
!
    dj = df(1, 1)*(df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2))-df(2, 1)*(df(1, 2)*df(3, 3)-df(1, 3)*df(&
         &3, 2))+df(3, 1)*(df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2))
!
    jp = jm*dj
!
!
!  CALCUL DE ETR
! ---------------
!
!    CALCUL DE BE EN T-
    do ij = 1, 6
        bem(ij) = kr(ij)-2*em(ij)
    end do
!
!
!    CALCUL PDF(AB,KL) = DF(A,K)*DF(B,L) SYMETRISE ET RACINE DE 2
    do ij = 1, 6
        i = ind1(ij)
        j = ind2(ij)
        do kl = 1, 6
            k = ind1(kl)
            l = ind2(kl)
            pdf(ij, kl) = rc(ij)*rc(kl)*(df(i, k)*df(j, l)+df(j, k)*df(i, l))/2
        end do
    end do
!
!
!    CALCUL DE BE TRIAL : BETR(AB) = PDF(AB,IJ):BEM(IJ)  ET  E TRIAL
    do ij = 1, 6
        b_n = to_blas_int(6)
        b_incx = to_blas_int(6)
        b_incy = to_blas_int(1)
        betr(ij) = ddot(b_n, pdf(ij, 1), b_incx, bem, b_incy)
        etr(ij) = (kr(ij)-betr(ij))/2
    end do
!
!
!    CALCUL DES INVARIANTS DE E TRIAL
    tretr = etr(1)+etr(2)+etr(3)
    do ij = 1, 6
        dvetr(ij) = etr(ij)-tretr/3.d0*kr(ij)
    end do
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    eqetr = sqrt(1.5d0*ddot(b_n, dvetr, b_incx, dvetr, b_incy))
!
end subroutine
