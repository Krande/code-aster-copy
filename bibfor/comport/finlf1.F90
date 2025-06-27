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

subroutine finlf1(ndim, delta, plasti, alpha, kn, kt, &
                  mu, Bn, Bt, m1, m2, cbar, D1, res)
!
    implicit none
!
    integer(kind=8) :: ndim
    real(kind=8) :: delta(ndim), plasti(ndim), alpha
    real(kind=8) :: kn, kt, mu, Bn, Bt, m1, m2, cbar, D1
    real(kind=8) :: res
!-----------------------------------------------------------------------
! IN : NDIM - NOMBRE DE COMPOSANTES (2 OU 3)
! IN : DELTA - SAUT DE DEPLACEMENT
! IN : PLASTI - VECTEUR DE DEPLACEMENT PLASTIQUE (INSTANT PRECEDENT)
! IN : ALPHA - VARIABLE D'ENDOMMAGEMENT (INSTANT PRESENT)
! IN : MU - COEFFICIENT DE FROTTEMENT
! IN : KN, KT - COEFFICIENTS DE RIGIDITE NORMALE ET TANGENTIELLE
! IN : RN, RT, M1, M2 - COEFFICIENTS DES FONCTIONS D'ENDOMMAGEMENT
! IN : D1 - ENERGIE MAXIMALE DISSIPE PAR ENDOMMAGEMENT
! OUT : RES - VALUE OF THE FUNCTION
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: wvec(ndim-1), wnorm, frwvec(ndim-1)
    real(kind=8) :: part1, part2, num1, den1, disc, disc2, den2
!
    wnorm = 0.d0
    do i = 1, ndim-1
        wvec(i) = kt*(delta(i+1)-plasti(i+1))*alpha**m2-Bt*plasti(i+1)*(1.d0-alpha)**m1
        wnorm = wnorm+wvec(i)**2
    end do
    wnorm = sqrt(wnorm)
    if (wnorm .le. 0.d0) then
        do i = 1, ndim-1
            frwvec(i) = 0.d0
        end do
    else
        do i = 1, ndim-1
            frwvec(i) = wvec(i)/wnorm
        end do
    end if
!
    num1 = wnorm+mu*kn*(delta(1)-plasti(1))*alpha**m2-mu*Bn*plasti(1)*(1-alpha)**m1-cbar*alpha**m2
    den1 = (mu**2*kn+kt)*alpha**m2+(mu**2*Bn+Bt)*(1.d0-alpha)**m1
    part1 = num1/den1*(1.d0-alpha)**((m1-1)/2.d0)
!
    disc = mu*Bn*plasti(1)
    do i = 1, ndim-1
        disc = disc+Bt*plasti(i+1)*frwvec(i)
    end do
    disc = disc*((m2-m1)*alpha-m2)
    disc = disc**2
    disc = disc*(1-alpha)**(m1-1)
    disc2 = Bn*plasti(1)**2
    do i = 2, ndim
        disc2 = disc2+Bt*plasti(i)**2
    end do
    disc2 = disc2*((m2-m1)*alpha-m2)*(1.d0-alpha)**(m1-1)+2.d0*D1*alpha**(m2+1)
    disc2 = disc2*(mu**2*Bn+Bt)*((m2-m1)*alpha-m2)
    disc = disc-disc2
!
    part2 = -mu*Bn*plasti(1)
    do i = 1, ndim-1
        part2 = part2-Bt*plasti(i+1)*frwvec(i)
    end do
    part2 = part2/(mu**2*Bn+Bt)*(1.d0-alpha)**((m1-1)/2.d0)
    den2 = (mu**2*Bn+Bt)*((m1-m2)*alpha+m2)
    part2 = part2+sqrt(disc)/den2
!
    res = (part1-part2)/sqrt(2.d0*D1/(mu**2*Bn+Bt)/m1)
!
end subroutine
