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
subroutine lcmhdd(necoul, necris, nbsys, nbcoef, coefh, &
                  nsg, hsr)
    implicit none
! person_in_charge: jean-michel.proix at edf.fr
!     ----------------------------------------------------------------
!     MONOCRISTAL : CALCUL DE LA MATRICE D'INTERACTION HSR POUR DD_CFC
!     ----------------------------------------------------------------
!     IN  NBSYS  :  NOMBRE DE SYSTEMES DE GLISSEMENT
!         NBCOEF  :  NOMBRE DE COEFFICIENTS
!         COEFH  :  COEFFICIENTS H1 A H6
!     OUT HSR    :  MATRICE D'INTERACTION
!     ----------------------------------------------------------------
#include "asterfort/lcicma.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbcoef, nbsys, i, j, nn(12), idbg, nsg
    real(kind=8) :: coefh(6), hsr(nsg, nsg), hgm(12, 12)
    real(kind=8) :: a1(3, 3), a2(3, 3), a3(3, 3), a4(3, 3), a0(3, 3), a5(3, 3)
    real(kind=8) :: a6(3, 3)
    real(kind=8) :: aetoil, acolin, agliss, alomer, ahirth, c0, c1, c2, c3, c4
    real(kind=8) :: c5
    character(len=16) :: necris, necoul
    data nn/7, 9, 8, 2, 1, 3, 12, 11, 10, 5, 4, 6/
!     ----------------------------------------------------------------
    idbg = 0
    if (necris(1:11) .eq. 'MONO_DD_CFC') then
!
!  MATRICE D INTERACTION (12*12): 5 COEFFICIENTS DD_CFC
!  DEFINITION SELON G.MONET
        if (nbsys .ne. 12) then
            call utmess('F', 'COMPOR1_24')
        end if
        aetoil = coefh(1)
        acolin = coefh(2)
        agliss = coefh(3)
        alomer = coefh(4)
        ahirth = coefh(5)
!
        call r8inir(3*3, aetoil, a0, 1)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 1)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 4)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 7)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 10)
!
        a1(1, 1) = acolin
        a1(1, 2) = agliss
        a1(1, 3) = agliss
        a1(2, 1) = agliss
        a1(2, 2) = ahirth
        a1(2, 3) = alomer
        a1(3, 1) = agliss
        a1(3, 2) = alomer
        a1(3, 3) = ahirth
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 1)
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 4)
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 7)
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 10)
!
        a2(1, 1) = ahirth
        a2(1, 2) = agliss
        a2(1, 3) = alomer
        a2(2, 1) = agliss
        a2(2, 2) = acolin
        a2(2, 3) = agliss
        a2(3, 1) = alomer
        a2(3, 2) = agliss
        a2(3, 3) = ahirth
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 1)
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 7)
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 4)
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 10)
!
        a3(1, 1) = ahirth
        a3(1, 2) = alomer
        a3(1, 3) = agliss
        a3(2, 1) = alomer
        a3(2, 2) = ahirth
        a3(2, 3) = agliss
        a3(3, 1) = agliss
        a3(3, 2) = agliss
        a3(3, 3) = acolin
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 4)
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 7)
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 1)
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 10)
!
        do i = 1, 12
            do j = 1, 12
                hsr(nn(i), nn(j)) = hgm(i, j)
            end do
        end do
!
    else if (necris(1:10) .eq. 'MONO_DD_CC') then
!
!
!  MATRICE D INTERACTION (12*12): 5 COEFFICIENTS DD_CFC
!  DEFINITION SELON G.MONET
        if (nbsys .ne. 12) then
            call utmess('F', 'COMPOR1_24')
        end if
        c0 = coefh(1)
        c1 = coefh(2)
        c2 = coefh(3)
        c3 = coefh(4)
        c4 = coefh(5)
        c5 = coefh(6)
!
        call r8inir(3*3, c1, a0, 1)
        do i = 1, 3
            a0(i, i) = c0
        end do
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 1)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 4)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 7)
        call lcicma(a0, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 10)
!
        a1(1, 1) = c4
        a1(1, 2) = c3
        a1(1, 3) = c2
        a1(2, 1) = c3
        a1(2, 2) = c5
        a1(2, 3) = c3
        a1(3, 1) = c2
        a1(3, 2) = c3
        a1(3, 3) = c4
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 7)
        call lcicma(a1, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 10)
!
        a2(1, 1) = c4
        a2(1, 2) = c2
        a2(1, 3) = c3
        a2(2, 1) = c2
        a2(2, 2) = c4
        a2(2, 3) = c3
        a2(3, 1) = c3
        a2(3, 2) = c3
        a2(3, 3) = c5
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 1)
        call lcicma(a2, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 10)
!
        a3(1, 1) = c5
        a3(1, 2) = c3
        a3(1, 3) = c3
        a3(2, 1) = c3
        a3(2, 2) = c4
        a3(2, 3) = c2
        a3(3, 1) = c3
        a3(3, 2) = c2
        a3(3, 3) = c4
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 1)
        call lcicma(a3, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 7)
!
!
        a4(1, 1) = c2
        a4(1, 2) = c3
        a4(1, 3) = c4
        a4(2, 1) = c3
        a4(2, 2) = c5
        a4(2, 3) = c3
        a4(3, 1) = c4
        a4(3, 2) = c3
        a4(3, 3) = c2
        call lcicma(a4, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    1, 4)
        call lcicma(a4, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 1)
!
        a5(1, 1) = c3
        a5(1, 2) = c3
        a5(1, 3) = c5
        a5(2, 1) = c2
        a5(2, 2) = c4
        a5(2, 3) = c3
        a5(3, 1) = c4
        a5(3, 2) = c2
        a5(3, 3) = c3
        call lcicma(a5, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    10, 4)
        call lcicma(a5, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 10)
!
        a6(1, 1) = c3
        a6(1, 2) = c2
        a6(1, 3) = c4
        a6(2, 1) = c3
        a6(2, 2) = c4
        a6(2, 3) = c2
        a6(3, 1) = c5
        a6(3, 2) = c3
        a6(3, 3) = c3
        call lcicma(a6, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    7, 4)
        call lcicma(a6, 3, 3, 3, 3, &
                    1, 1, hgm, 12, 12, &
                    4, 7)
!
        do i = 1, 12
            do j = 1, 12
                hsr(i, j) = hgm(i, j)
            end do
        end do
        if (idbg .eq. 1) then
            write (6, *) 'MATRICE D INTERACTION POUR', necris
            do i = 1, 12
                write (6, '(12(1X,E11.4))') (hgm(i, j), j=1, 12)
            end do
        end if
!
    end if
!
!
end subroutine
