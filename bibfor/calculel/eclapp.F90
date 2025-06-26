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
subroutine eclapp(ndim, nno2, lonmin, coor)
!
    implicit none
#include "asterc/matfpe.h"
#include "asterfort/r8inir.h"
#include "blas/dnrm2.h"
    integer(kind=8) :: ndim, nno2
    real(kind=8) :: coor(ndim, nno2), lonmin
!
! ----------------------------------------------------------------------
!        CONTROLE DE L'APPLATISSEMENT DES ELEMENTS AVEC ECLA_PG
! ----------------------------------------------------------------------
! IN  NDIM    I  DIMENSION DE L'ESPACE
! IN  NNO2    I  NOMBRE DE NOEUDS DE L'ELEMENT
! IN  LONMIN  R  LONGUEUR MINIMALE (LONGUEUR DES MEDIANES)
! VAR COOR    R  COORDONNEES DES NOEUDS
!                 IN  -> VALEURS INITIALES
!                 OUT -> VALEURS APRES CORRECTION (AFFINITE)
! ----------------------------------------------------------------------
!
    real(kind=8) :: corr, l1, l2, d1(3), d2(3), prec
    integer(kind=8) :: i
    blas_int :: b_incx, b_n
    parameter(prec=1.d-5)
! ----------------------------------------------------------------------
    call matfpe(-1)
!
    if (nno2 .eq. 4) then
!
        do i = 1, ndim
            d1(i) = (coor(i, 2)+coor(i, 3)-coor(i, 1)-coor(i, 4))/2
            d2(i) = (coor(i, 3)+coor(i, 4)-coor(i, 1)-coor(i, 2))/2
        end do
!
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        l1 = dnrm2(b_n, d1, b_incx)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        l2 = dnrm2(b_n, d2, b_incx)
!
!      ELEMENTS PLATS
        if (min(l1, l2)/max(l1, l2) .lt. prec) then
            if (l1 .lt. l2) then
                call r8inir(ndim, 0.d0, d1, 1)
                d1(1) = d2(2)
                d1(2) = -d2(1)
                b_n = to_blas_int(ndim)
                b_incx = to_blas_int(1)
                corr = lonmin/2.d0/dnrm2(b_n, d1, b_incx)
                do i = 1, ndim
                    coor(i, 1) = coor(i, 1)-corr*d1(i)
                    coor(i, 2) = coor(i, 2)+corr*d1(i)
                    coor(i, 3) = coor(i, 3)+corr*d1(i)
                    coor(i, 4) = coor(i, 4)-corr*d1(i)
                end do
            else
                call r8inir(ndim, 0.d0, d2, 1)
                d2(1) = -d1(2)
                d2(2) = d1(1)
                b_n = to_blas_int(ndim)
                b_incx = to_blas_int(1)
                corr = lonmin/2.d0/dnrm2(b_n, d2, b_incx)
                do i = 1, ndim
                    coor(i, 1) = coor(i, 1)-corr*d2(i)
                    coor(i, 2) = coor(i, 2)-corr*d2(i)
                    coor(i, 3) = coor(i, 3)+corr*d2(i)
                    coor(i, 4) = coor(i, 4)+corr*d2(i)
                end do
            end if
            goto 999
        end if
!
!      ELEMENTS EPAIS
        if (l1 .lt. lonmin) then
            corr = lonmin/l1
            do i = 1, ndim
                coor(i, 1) = corr*coor(i, 1)+(1-corr)/2*(coor(i, 1)+coor(i, 2))
                coor(i, 2) = corr*coor(i, 2)+(1-corr)/2*(coor(i, 1)+coor(i, 2))
                coor(i, 3) = corr*coor(i, 3)+(1-corr)/2*(coor(i, 3)+coor(i, 4))
                coor(i, 4) = corr*coor(i, 4)+(1-corr)/2*(coor(i, 3)+coor(i, 4))
            end do
        end if
!
        if (l2 .lt. lonmin) then
            corr = lonmin/l2
            do i = 1, ndim
                coor(i, 1) = corr*coor(i, 1)+(1-corr)/2*(coor(i, 1)+coor(i, 4))
                coor(i, 2) = corr*coor(i, 2)+(1-corr)/2*(coor(i, 2)+coor(i, 3))
                coor(i, 3) = corr*coor(i, 3)+(1-corr)/2*(coor(i, 2)+coor(i, 3))
                coor(i, 4) = corr*coor(i, 4)+(1-corr)/2*(coor(i, 1)+coor(i, 4))
            end do
        end if
!
    end if
!
999 continue
!
    call matfpe(1)
!
end subroutine
