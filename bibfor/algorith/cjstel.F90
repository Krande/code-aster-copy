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
subroutine cjstel(mod, mater, sig, hook)
    implicit none
#include "asterfort/utmess.h"
!     CALCUL DE LA MATRICE DE RIGIDITE ELASTIQUE DE LA LOI CJS
!     IN   MOD     :  MODELISATION
!          MATER   :  COEFFICIENTS MATERIAU
!          SIG     :  CONTRAINTES
!     OUT  HOOK    :  OPERATEUR RIGIDITE ELASTIQUE
!       ----------------------------------------------------------------
!
    integer(kind=8) :: ndt, ndi
!
    real(kind=8) :: sig(6), hook(6, 6), mater(14, 2), i1, coef
    real(kind=8) :: e, nu, al, la, mu
    real(kind=8) :: un, d12, zero, deux, trois, qinit
    integer(kind=8) :: i, j
!
    character(len=8) :: mod
!
    common/tdim/ndt, ndi
!
    data d12/.5d0/
    data un/1.d0/
    data zero/0.d0/
    data deux/2.d0/
    data trois/3.d0/
!
!       ----------------------------------------------------------------
!
!--->   CALCUL PREMIER INVARIANT DES CONTRAINTES
    qinit = mater(13, 2)
    i1 = zero
    do i = 1, ndi
        i1 = i1+sig(i)
    end do
!
!
!--->   CALCUL DES COEF. UTILES
    if ((i1+qinit) .eq. zero) then
        coef = un
    else
        coef = ((i1+qinit)/trois/mater(12, 2))**mater(3, 2)
    end if
    e = mater(1, 1)*coef
    nu = mater(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
!
!
!--->   OPERATEUR DE RIGIDITE
!
!
! - 3D/DP/AX
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        do i = 1, ndi
            do j = 1, ndi
                if (i .eq. j) hook(i, j) = al
                if (i .ne. j) hook(i, j) = la
            end do
        end do
        do i = ndi+1, ndt
            do j = ndi+1, ndt
                if (i .eq. j) hook(i, j) = deux*mu
            end do
        end do
!
! - CP/1D
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
    end if
!
end subroutine
