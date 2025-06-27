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

subroutine finlf2(ndim, delta, alpha, kn, kt, &
                  mu, Bn, Bt, m1, m2, cbar, d1, res)
!
    implicit none
!
    integer(kind=8) :: ndim
    real(kind=8) :: delta(ndim), alpha
    real(kind=8) :: kn, kt, mu, Bn, Bt, m1, m2, cbar, d1
    real(kind=8) :: res
!-----------------------------------------------------------------------
! IN : NDIM - NOMBRE DE COMPOSANTES (2 OU 3)
! IN : DELTA - SAUT DE DEPLACEMENT
! IN : ALPHA - VARIABLE D'ENDOMMAGEMENT (INSTANT PRESENT)
! IN : MU - COEFFICIENT DE FROTTEMENT
! IN : KN, KT - COEFFICIENTS DE RIGIDITE NORMALE ET TANGENTIELLE
! IN : RN, RT, M1, M2 - COEFFICIENTS DES FONCTIONS D'ENDOMMAGEMENT
! IN : D1 - ENERGIE MAXIMALE DISSIPE PAR ENDOMMAGEMENT
! OUT : RES - VALUE OF THE FUNCTION
!-----------------------------------------------------------------------
    integer(kind=8) :: i
!
    res = Bn*((mu*kn*delta(1)-cbar)/(mu*(kn*alpha**m2+Bn*(1-alpha)**m1)))**2
    do i = 2, ndim
        res = res+Bt*(kt*delta(i)/(kt*alpha**m2+Bt*(1-alpha)**m1))**2
    end do
    res = res*((m1-m2)*alpha+m2)*(1-alpha)**(m1-1)
    res = res/d1-2*(alpha**(1-m2))
end subroutine
