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
subroutine ceobfb(bm, epsm, lambda, mu, ecrob, &
                  bdim, fb, nofbm, fbm)
! person_in_charge: ludovic.idoux at edf.fr
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/diago2.h"
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    real(kind=8) :: epsm(6), bm(6), fb(6), fbm(6), nofbm
    real(kind=8) :: lambda, mu, ecrob
    integer(kind=8) :: bdim
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DU MODELE D'ENDOMMAGEMENT ANISOTROPE
!     ROUTINE DE CALCUL DE LA FORCE THERMODYNAMIQUE FB
!
!  IN BM     : TENSEUR D'ENDOMMAGEMENT DE TRACTION
!  IN EPSM   : TENSEUR DE DEFORMATION
!  IN LAMBDA : /
!  IN MU     : / COEFFICIENTS DE LAME
!  IN ECROB  : PARAMETRE DU MODELE
!  IN BDIM   : DIMENSION DE L ESPACE
!
! OUT FB     : FORCE THERMODYNAMIQUE
! OUT FBM    : PARTIE POSITIVE DE FB
! OUT NOFBM  : NORME DE FBM
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k
    integer(kind=8) :: t(3, 3), r(2, 2)
!
    real(kind=8) :: cc(6), cpe(6), ccp(6), eps(6), b(6), fbs(3)
    real(kind=8) :: deux, treb, kron(6)
    real(kind=8) :: vecc(3, 3), valcc(3), vecfb(3, 3), valfb(3)
    real(kind=8) :: vecfbs(2, 2), valfbs(2)
!
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
    deux = 2.d0
!
    do i = 1, 6
        b(i) = bm(i)
        eps(i) = epsm(i)
    end do
!
! CALCUL DE FB
!
    call r8inir(6, 0.d0, cc, 1)
!
    do i = 1, 3
        do j = i, 3
            do k = 1, 3
                cc(t(i, j)) = cc(t(i, j))+b(t(i, k))*eps(t(k, j))+b(t(j, k)) &
                              *eps(t(k, i))
            end do
        end do
    end do
    call diago3(cc, vecc, valcc)
    call r8inir(6, 0.d0, ccp, 1)
    call r8inir(6, 0.d0, cpe, 1)
    do i = 1, 3
        if (valcc(i) .lt. 0.d0) then
            valcc(i) = 0.d0
        end if
    end do
    do i = 1, 3
        do j = i, 3
            do k = 1, 3
                ccp(t(i, j)) = ccp(t(i, j))+vecc(i, k)*valcc(k)*vecc(j, k)
            end do
        end do
    end do
    do i = 1, 3
        do j = i, 3
            do k = 1, 3
                cpe(t(i, j)) = cpe(t(i, j))+ccp(t(i, k))*eps(t(k, j))+ &
                               ccp(t(j, k))*eps(t(k, i))
            end do
        end do
    end do
!
    call r8inir(6, 0.d0, fb, 1)
    treb = 0.d0
    do i = 1, 3
        treb = treb+cc(i)/deux
    end do
    if (treb .gt. 0.d0) then
        do i = 1, 6
            fb(i) = -lambda*treb*eps(i)
        end do
    end if
    do i = 1, 6
        fb(i) = fb(i)-mu/deux*cpe(i)+ecrob*(kron(i)-b(i))
    end do
!
! CALCUL DE LA PARTIE POSITIVE DE FBM ET DE SA NORME NOFB
!
    call r8inir(6, 0.d0, fbm, 1)
    if (bdim .eq. 3) then
        call diago3(fb, vecfb, valfb)
        nofbm = 0.d0
!
        do i = 1, 3
            if (valfb(i) .gt. 0.d0) then
                valfb(i) = 0.d0
            end if
            nofbm = nofbm+valfb(i)*valfb(i)
        end do
!
        do i = 1, 3
            do j = i, 3
                do k = 1, 3
                    fbm(t(i, j)) = fbm(t(i, j))+vecfb(i, k)*valfb(k)*vecfb( &
                                   j, k)
                end do
            end do
        end do
!
    else if (bdim .eq. 2) then
        r(1, 1) = 1
        r(2, 2) = 2
        r(1, 2) = 3
        r(2, 1) = 3
        fbs(1) = fb(1)
        fbs(2) = fb(2)
        fbs(3) = fb(4)
!
        call diago2(fbs, vecfbs, valfbs)
!
        nofbm = 0.d0
        do i = 1, 2
            if (valfbs(i) .gt. 0.d0) then
                valfbs(i) = 0.d0
            end if
            nofbm = nofbm+valfbs(i)*valfbs(i)
        end do
!
        do i = 1, 2
            do j = i, 2
                do k = 1, 2
                    fbm(r(i, j)) = fbm(r(i, j))+vecfbs(i, k)*valfbs(k)* &
                                   vecfbs(j, k)
                end do
            end do
        end do
!
    else if (bdim .eq. 1) then
        if (fb(1) .lt. 0.d0) then
            fbm(1) = fb(1)
        end if
        nofbm = fbm(1)**2
!
    end if
!
    if (abs(nofbm) .lt. r8prem()) then
        nofbm = 0.d0
    end if
!
end subroutine
