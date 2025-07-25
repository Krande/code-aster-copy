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
subroutine gdsmci(fm, df, em)
!
!
    implicit none
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    real(kind=8) :: fm(3, 3), df(3, 3), em(6)
!
! ----------------------------------------------------------------------
!       INTEGRATION DES LOIS EN GRANDES DEFORMATIONS SIMO-MIEHE
!                  CALCUL DES ELEMENTS CINEMATIQUES
! ----------------------------------------------------------------------
! IN  FM    : DEFORMATION AU DEBUT DU PAS DE TEMPS
! IN  DF    : INCREMENT DE DEFORMATION PENDANT LE PAS DE TEMPS
! IN  EM    : DEFORMATION ELASTIQUE AU DEBUT DU PAS DE TEMPS (VIM)
! ----------------------------------------------------------------------
! COMMON GRANDES DEFORMATIONS SIMO - MIEHE
!
    integer(kind=8) :: ind(3, 3), ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6), id(6, 6)
    real(kind=8) :: bem(6), betr(6), dvbetr(6), eqbetr, trbetr
    real(kind=8) :: jp, dj, jm, dfb(3, 3)
    real(kind=8) :: djdf(3, 3), dbtrdf(6, 3, 3)
!
    common/gdsmc/&
     &            bem, betr, dvbetr, eqbetr, trbetr,&
     &            jp, dj, jm, dfb,&
     &            djdf, dbtrdf,&
     &            kr, id, rac2, rc, ind, ind1, ind2
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ij, kl, i, j, k, l
    real(kind=8) :: pdf(6, 6)
!
    real(kind=8) :: e(3, 3), dete, fp(3, 3)
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
!
!
! 1 - CALCUL DES JACOBIENS ET DE DF-BARRE
! ----------------------------------------
!
    jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2))-fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(3&
         &, 2))+fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
!
    dj = df(1, 1)*(df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2))-df(2, 1)*(df(1, 2)*df(3, 3)-df(1, 3)*df(3&
         &, 2))+df(3, 1)*(df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2))
!
    jp = jm*dj
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, df, b_incx, dfb, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    call dscal(b_n, dj**(-1.d0/3.d0), dfb, b_incx)
!
!
!
!  CALCUL DE BETR
! ---------------
!
!    CALCUL DE BE-BARRE EN T-
    do ij = 1, 6
        bem(ij) = (kr(ij)-2*em(ij))/jm**(2.d0/3.d0)
    end do
!
!
!    CALCUL PDF(IJ,KL) = DFB(I,K)*DFB(J,L) SYMETRISE ET RACINE DE 2
    do ij = 1, 6
        i = ind1(ij)
        j = ind2(ij)
        do kl = 1, 6
            k = ind1(kl)
            l = ind2(kl)
            pdf(ij, kl) = rc(ij)*rc(kl)*(dfb(i, k)*dfb(j, l)+dfb(j, k)*dfb(i, l))/2.d0
        end do
    end do
!
!
!    CALCUL DE BE TRIAL : BETR(AB) = PDF(AB,IJ):BEM(IJ)
    do ij = 1, 6
        b_n = to_blas_int(6)
        b_incx = to_blas_int(6)
        b_incy = to_blas_int(1)
        betr(ij) = ddot(b_n, pdf(ij, 1), b_incx, bem, b_incy)
    end do
!
!
    do i = 1, 3
        do j = 1, 3
            fp(i, j) = 0
            do k = 1, 3
                fp(i, j) = fp(i, j)+df(i, k)*fm(k, j)
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            e(i, j) = 0
            do k = 1, 3
                e(i, j) = e(i, j)+fp(i, k)*fp(j, k)
            end do
        end do
    end do
!
    dete = e(1, 1)*(e(2, 2)*e(3, 3)-e(2, 3)*e(3, 2))-e(2, 1)*(e(1, 2)*e(3, 3)-e(1, 3)*e(3, 2))+e(&
           &3, 1)*(e(1, 2)*e(2, 3)-e(1, 3)*e(2, 2))
!
    do i = 1, 3
        do j = 1, 3
            e(i, j) = e(i, j)/dete**(1.d0/3.d0)
        end do
    end do
!
!
! 3.4 - CALCUL DE LA PARTIE DEVIATORIQUE ET INVARIANTS DE BETR EN T+
!
    trbetr = betr(1)+betr(2)+betr(3)
    do ij = 1, 6
        dvbetr(ij) = betr(ij)-trbetr/3.d0*kr(ij)
    end do
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    eqbetr = sqrt(1.5d0*ddot(b_n, dvbetr, b_incx, dvbetr, b_incy))
!
end subroutine
