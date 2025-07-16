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
subroutine gdcltg(df, e)
!
!
    implicit none
#include "asterfort/r8inir.h"
#include "blas/dscal.h"
    real(kind=8) :: df(3, 3), e(6)
!
! ----------------------------------------------------------------------
!       INTEGRATION DES LOIS EN GRANDES DEFORMATIONS CANO-LORENTZ
!   CALCUL DES DERIVEES PAR RAPPORT A UNE VARIATION DE LA DEFORMATION
! ----------------------------------------------------------------------
! IN  DF      INCREMENT DE DEFORMATION PENDANT LE PAS DE TEMPS
! IN  E       DEFORMATION ELASTIQUE (XX,YY,ZZ,RAC2*XY,RAC2*XZ,RAC2*YZ)
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
!
    integer(kind=8) :: ij, kl, i, j, l, il, jl, ind12(3, 3)
    real(kind=8) :: tre, coef, betr(6)
    blas_int :: b_incx, b_n
    data ind12/1, 4, 5, 4, 2, 6, 5, 6, 3/
! ----------------------------------------------------------------------
!
!
!
!  CALCUL DE LA DERIVEE DES JACOBIENS DJ / DF = DJDF
! ----------------------------------------------------
!
    call r8inir(9, 0.d0, djdf, 1)
    djdf(1, 1) = jp
    djdf(2, 2) = jp
    djdf(3, 3) = jp
!
!
!  CALCUL DE LA DERIVEE DE DETR / DF : DETRDF(AB,P,Q)
! ----------------------------------------------------
!
    do ij = 1, 6
        betr(ij) = kr(ij)-2*etr(ij)
    end do
    call r8inir(54, 0.d0, detrdf, 1)
!
    do ij = 1, 6
        i = ind1(ij)
        j = ind2(ij)
        do l = 1, 3
            il = ind12(i, l)
            jl = ind12(j, l)
            detrdf(ij, i, l) = detrdf(ij, i, l)-0.5d0*rc(ij)*betr(jl)/rc(jl)
            detrdf(ij, j, l) = detrdf(ij, j, l)-0.5d0*rc(ij)*betr(il)/rc(il)
        end do
    end do
!
!
!  DERIVEE PARTIELLE DE TAU PAR RAPPORT A E
! --------------------------------------------
!
    call r8inir(36, 0.d0, dtaude, 1)
!
!    TERME D(E.E)/DE  DE DTAU/DE
    dtaude(1, 1) = 2*e(1)
    dtaude(2, 2) = 2*e(2)
    dtaude(3, 3) = 2*e(3)
    dtaude(4, 1) = e(4)
    dtaude(4, 2) = e(4)
    dtaude(5, 1) = e(5)
    dtaude(5, 3) = e(5)
    dtaude(6, 2) = e(6)
    dtaude(6, 3) = e(6)
    dtaude(4, 4) = e(1)+e(2)
    dtaude(5, 5) = e(1)+e(3)
    dtaude(6, 6) = e(2)+e(3)
    dtaude(5, 4) = e(6)/rac2
    dtaude(6, 4) = e(5)/rac2
    dtaude(6, 5) = e(4)/rac2
!
    do ij = 1, 6
        do kl = ij+1, 6
            dtaude(ij, kl) = dtaude(kl, ij)
        end do
    end do
    b_n = to_blas_int(36)
    b_incx = to_blas_int(1)
    call dscal(b_n, 2*deuxmu, dtaude, b_incx)
!
!    TERME EN (2E-1) X 1  DE DTAU / DE
    do ij = 1, 6
        do kl = 1, 3
            dtaude(ij, kl) = dtaude(ij, kl)+lambda*(2*e(ij)-kr(ij))
        end do
    end do
!
!    TERME EN ID  DE DTAU/DE
    tre = e(1)+e(2)+e(3)
    coef = 2*(lambda*tre+cother)-deuxmu
    dtaude(1, 1) = dtaude(1, 1)+coef
    dtaude(2, 2) = dtaude(2, 2)+coef
    dtaude(3, 3) = dtaude(3, 3)+coef
    dtaude(4, 4) = dtaude(4, 4)+coef
    dtaude(5, 5) = dtaude(5, 5)+coef
    dtaude(6, 6) = dtaude(6, 6)+coef
!
end subroutine
