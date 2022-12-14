! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine brksec(h66, bt3, bc, nu, e,&
                  s3)
!
!      H66 MATRICE SECANTE
!      B  INDICE DE FISSURATION  EN BASE PRINCIPALE D ENDO
!      BC INDICE DE FISSURATION COMP
!      E,NU ELAST MATR SAIN
!
!     CALCUL DE LA MATRICE SECANTE ORTHOTROPE EN
!       BASE PRINCIPALE D ENDOMMAGEMENT
    implicit none
#include "asterc/r8prem.h"
    real(kind=8) :: h66(6, 6), b(3), s3(3), nu, e, bt3(3)
    integer :: i, j
    real(kind=8) :: bc, t1, t17, t2, t21, t24, t26
    real(kind=8) :: t31, t5, t7, t8
!-----------------------------------------------------------------------
!
!     MISE  AZERO
    do i = 1, 6
        do j = 1, 6
            h66(i,j)=0.d0
        end do
    end do
!
!
!     PRISE EN COMPTE DU CARACTERE UNILATERAL
    do i = 1, 3
        if (s3(i) .le. r8prem()*e) then
            bt3(i)=0.d0
        else
            if (bt3(i) .gt. 2.3d0) then
                bt3(i)=2.3d0
            endif
        endif
        b(i)=exp(bt3(i))
    end do
!
!     CARRé SUPERIEUR RELIANT LES CONTRAINTES NORMALES DANS
!       LA MATRICE SECANTE
    t1 = b(2)
    t2 = b(3)
    t5 = nu ** 2
    t7 = b(1)
    t8 = t7 * t1
    t17 = 0.1d1 / ( -0.1d1 * t8 * t2 + t7 * t5 + t5 * t2 + 0.2d1 * t5 * nu + t5 * t1)
    t21 = (t2 + nu) * nu * t17
    t24 = (nu + t1) * nu * t17
    t26 = 0.1d1 * t5
    t31 = (t7 + nu) * nu * t17
    h66(1,1) = (-0.1d1 * t1 * t2 + t5) * t17
    h66(1,2) = -t21
    h66(1,3) = -t24
    h66(2,1) = -t21
    h66(2,2) = -(t7 * t2 - t26) * t17
    h66(2,3) = -t31
    h66(3,1) = -t24
    h66(3,2) = -t31
    h66(3,3) = -(t8 - t26) * t17
!     SUITE DE LA MATRICE SECANTE (TERMES DE CISAILLEMENT)
    h66(4,4) = exp(-(bt3(1)+bt3(2)))/(1.d0+nu)
    h66(5,5) = exp(-(bt3(1)+bt3(3)))/(1.d0+nu)
    h66(6,6) = exp(-(bt3(2)+bt3(3)))/(1.d0+nu)
!
    do i = 1, 6
        do j = 1, 6
            h66(i,j)=e*h66(i,j)*exp(-bc)
        end do
    end do
end subroutine
