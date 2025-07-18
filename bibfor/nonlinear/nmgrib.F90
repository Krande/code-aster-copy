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
subroutine nmgrib(nno, geom, dff, dir11, lexc, &
                  vecn, b, jac, p)
! CALCUL DE LA MATRICE B ET JACOBIEN POUR LES GRILLES SECONDE GENERATION
! ----------------------------------------------------------------------
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/r8inir.h"
#include "asterfort/subaco.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
#include "asterfort/utmess.h"
    aster_logical :: lexc
    integer(kind=8) :: nno
    real(kind=8) :: geom(3, nno), dff(2, nno), dir11(3)
    real(kind=8) :: b(6, nno), vecn(3), p(3, 6)
    integer(kind=8) :: i, j, n, alpha, beta, gamma
    real(kind=8) :: cova(3, 3), metr(2, 2), jac, cnva(3, 2), a(2, 2), r1(3)
    real(kind=8) :: projn
    real(kind=8) :: mtemp(3, nno), denomi
!
    call subaco(nno, dff, geom, cova)
    call sumetr(cova, metr, jac)
    call subacv(cova, metr, jac, cnva, a)
!
    call r8inir(3, 0.d0, r1, 1)
    call r8inir(6*nno, 0.d0, b, 1)
!
    projn = 0.d0
!
    do j = 1, 3
        do i = 1, 2
            r1(i) = r1(i)+cova(j, i)*dir11(j)
        end do
        projn = projn+cova(j, 3)*dir11(j)
    end do
!
    denomi = (1.d0-projn*projn)
    if (abs(denomi) .le. r8prem()) then
        call utmess('F', 'ELEMENTS_3')
    end if
!
    do i = 1, 3
        do n = 1, nno
            do alpha = 1, 2
                do beta = 1, 2
                    do gamma = 1, 2
                       b(i, n) = b(i, n)+r1(alpha)*r1(gamma)*a(beta, gamma)*dff(beta, n)*cnva(i, al&
                                   &pha)/denomi
                    end do
                end do
            end do
        end do
    end do
!
    if (lexc) then
        do n = 1, nno
            do i = 1, 3
                mtemp(i, n) = b(i, n)
            end do
        end do
!
        call r8inir(18, 0.d0, p, 1)
        call r8inir(6*nno, 0.d0, b, 1)
!
        do i = 1, 3
            p(i, i) = 1.d0
        end do
        p(1, 5) = vecn(3)
        p(1, 6) = -vecn(2)
        p(2, 4) = -vecn(3)
        p(2, 6) = vecn(1)
        p(3, 4) = vecn(2)
        p(3, 5) = -vecn(1)
!
        do n = 1, nno
            do i = 1, 6
                do j = 1, 3
                    b(i, n) = b(i, n)+mtemp(j, n)*p(j, i)
                end do
            end do
        end do
!
    end if
end subroutine
