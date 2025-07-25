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
subroutine nuainr(method, np1, nx1, nc1, ic1, &
                  nuax1, nual1, nuav1, x2, dref, &
                  val2)
    implicit none
#include "asterf_types.h"
#include "asterfort/mgauss.h"
#include "asterfort/utmess.h"
    character(len=*) :: method
    integer(kind=8) :: nx1, np1, ic1, nc1
    real(kind=8) :: nuax1(*), nuav1(*), x2(nx1), dref, dref2, val2
    aster_logical :: nual1(*)
!
!  BUT : INTERPOLER LA VALEUR VAL2 DU NUAGE NUAG1 SUR LE POINT DE
!        COORDONNEES X2
!
! IN  METHOD   : METHODE D'INTERPOLATION: 'NUAGE_DEG_0' OU 'NUAGE_DEG_1'
! IN  NP1      : NOMBRE DE POINTS DU NUAGE NUAG1
! IN  NX1      : NOMBRE DE COORDONNEES DES POINTS DE NUAG1 (ET X2)
! IN  NC1      : NOMBRE DE CMPS DES POINTS DE NUAG1
! IN  IC1      : NUMERO DE LA CMP A INTERPOLER DANS NUAG1
! IN  NUAX1    : OBJET .NUAX  DE NUAG1
! IN  NUAL1    : OBJET .NUAL  DE NUAG1
! IN  NUAV1    : OBJET .NUAV  DE NUAG1
! IN  X2       : COORDONNEES DU POINT QU L'ON CHERCHE A INTERPOLER
! IN  DREF     : DISTANCE DE REFERENCE POUR LA METHODE D'INTERPOLATION
! OU  VAL2     : VALEUR INTERPOLEE
!
! VARIABLES LOCALES :
    integer(kind=8) :: ip1, ix1, i, j, iret
    real(kind=8) :: k0(1, 1), k1(2, 2), k2(3, 3), k3(4, 4), f(4)
    real(kind=8) :: w, x1, y1, z1, v1, d, det
    character(len=16) :: meth2
!
! PARAMETRES DE L'EXPONENTIELLE DONNANT LE POIDS DES POINTS :
    real(kind=8) :: alpha, beta
    data alpha, beta/.2d0, 0.75d0/
!
! -DEB
!
!     -- MISE A ZERO DE K ET F :
!     -------------------------
    if (nx1 .eq. 1) then
        do i = 1, nx1+1
            do j = 1, nx1+1
                k1(i, j) = 0.d0
            end do
        end do
    else if (nx1 .eq. 2) then
        do i = 1, nx1+1
            do j = 1, nx1+1
                k2(i, j) = 0.d0
            end do
        end do
    else if (nx1 .eq. 3) then
        do i = 1, nx1+1
            do j = 1, nx1+1
                k3(i, j) = 0.d0
            end do
        end do
    end if
!
    do i = 1, nx1+1
        f(i) = 0.d0
    end do
!
    dref2 = dref*alpha
!
    if (method .eq. 'NUAGE_DEG_0') then
!     ------------------------
        k0(1, 1) = 0.d0
        f(1) = 0.d0
!
        if (nx1 .eq. 1) then
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 1
                d = (nuax1((ip1-1)*1+1)-x2(1))**2
                w = exp(-(d/dref2)**beta)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k0(1, 1) = k0(1, 1)+w
                f(1) = f(1)+w*v1
1               continue
            end do
!
        else if (nx1 .eq. 2) then
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 2
                d = (nuax1((ip1-1)*2+1)-x2(1))**2+(nuax1((ip1-1)*2+2)- &
                                                   x2(2))**2
                w = exp(-(d/dref2)**beta)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k0(1, 1) = k0(1, 1)+w
                f(1) = f(1)+w*v1
2               continue
            end do
!
        else if (nx1 .eq. 3) then
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 3
                d = (nuax1((ip1-1)*3+1)-x2(1))**2+(nuax1((ip1-1)*3+2)- &
                                                   x2(2))**2+(nuax1((ip1-1)*3+3)-x2(3))**2
                w = exp(-(d/dref2)**beta)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k0(1, 1) = k0(1, 1)+w
                f(1) = f(1)+w*v1
3               continue
            end do
        end if
!
!
        val2 = f(1)/k0(1, 1)
!
!
    else if (method .eq. 'NUAGE_DEG_1') then
!     -----------------------------
        if (nx1 .eq. 1) then
!       --------------
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 11
                d = (nuax1((ip1-1)*1+1)-x2(1))**2
                w = exp(-(d/dref2)**beta)
!
                x1 = nuax1((ip1-1)*1+1)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k1(1, 1) = k1(1, 1)+w
                k1(1, 2) = k1(1, 2)+w*x1
                k1(2, 2) = k1(2, 2)+w*x1*x1
                f(1) = f(1)+w*v1
                f(2) = f(2)+w*v1*x1
11              continue
            end do
            k1(2, 1) = k1(1, 2)
!
        else if (nx1 .eq. 2) then
!       ------------------
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 21
                d = (nuax1((ip1-1)*2+1)-x2(1))**2+(nuax1((ip1-1)*2+2)- &
                                                   x2(2))**2
                w = exp(-(d/dref2)**beta)
!
                x1 = nuax1((ip1-1)*2+1)
                y1 = nuax1((ip1-1)*2+2)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k2(1, 1) = k2(1, 1)+w
                k2(1, 2) = k2(1, 2)+w*x1
                k2(1, 3) = k2(1, 3)+w*y1
                k2(2, 2) = k2(2, 2)+w*x1*x1
                k2(2, 3) = k2(2, 3)+w*x1*y1
                k2(3, 3) = k2(3, 3)+w*y1*y1
                f(1) = f(1)+w*v1
                f(2) = f(2)+w*v1*x1
                f(3) = f(3)+w*v1*y1
21              continue
            end do
            k2(2, 1) = k2(1, 2)
            k2(3, 1) = k2(1, 3)
            k2(3, 2) = k2(2, 3)
!
        else if (nx1 .eq. 3) then
!       ------------------
            do ip1 = 1, np1
                if (.not. nual1((ip1-1)*nc1+ic1)) goto 31
                d = (nuax1((ip1-1)*3+1)-x2(1))**2+(nuax1((ip1-1)*3+2)- &
                                                   x2(2))**2+(nuax1((ip1-1)*3+3)-x2(3))**2
                w = exp(-(d/dref2)**beta)
!
                x1 = nuax1((ip1-1)*3+1)
                y1 = nuax1((ip1-1)*3+2)
                z1 = nuax1((ip1-1)*3+3)
                v1 = nuav1((ip1-1)*nc1+ic1)
                k3(1, 1) = k3(1, 1)+w
                k3(1, 2) = k3(1, 2)+w*x1
                k3(1, 3) = k3(1, 3)+w*y1
                k3(1, 4) = k3(1, 4)+w*z1
                k3(2, 2) = k3(2, 2)+w*x1*x1
                k3(2, 3) = k3(2, 3)+w*x1*y1
                k3(2, 4) = k3(2, 4)+w*x1*z1
                k3(3, 3) = k3(3, 3)+w*y1*y1
                k3(3, 4) = k3(3, 4)+w*y1*z1
                k3(4, 4) = k3(4, 4)+w*z1*z1
                f(1) = f(1)+w*v1
                f(2) = f(2)+w*v1*x1
                f(3) = f(3)+w*v1*y1
                f(4) = f(4)+w*v1*z1
31              continue
            end do
            k3(2, 1) = k3(1, 2)
            k3(3, 1) = k3(1, 3)
            k3(3, 2) = k3(2, 3)
            k3(4, 1) = k3(1, 4)
            k3(4, 2) = k3(2, 4)
            k3(4, 3) = k3(3, 4)
        end if
        if (nx1 .eq. 1) then
            call mgauss('NFVP', k1, f, 2, 2, &
                        1, det, iret)
        else if (nx1 .eq. 2) then
            call mgauss('NFVP', k2, f, 3, 3, &
                        1, det, iret)
        else if (nx1 .eq. 3) then
            call mgauss('NFVP', k3, f, 4, 4, &
                        1, det, iret)
        end if
!
        val2 = f(1)
        do ix1 = 1, nx1
            val2 = val2+f(ix1+1)*x2(ix1)
        end do
!
    else
        meth2 = method
        call utmess('F', 'UTILITAI2_60', sk=meth2)
    end if
!
!
end subroutine
