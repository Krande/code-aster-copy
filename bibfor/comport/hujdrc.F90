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

subroutine hujdrc(k, mater, sig, vin, pst)
    implicit none
!   CALCUL DE PRODUIT SCALAIRE ENTRE LA NORME DES MECANISMES CYCLIQUES
!   DEVIATOIRES ET DES VECTEURS DE RAYONS ISSUS VARIABLES D'HISTOIRE
!   ET LA POSITION ACTUELLE DANS LE PLAN DEVIATOIRE K.
!
!   IN  K      :  PLAN DE PROJECTION (K = 1 A 3)
!       MATER  :  COEFFICIENTS MATERIAU A T+DT
!       VIN    :  VARIABLES INTERNES  A T
!       SIG    :  CONTRAINTE A T+DT
!
!   OUT PST    : PRODUIT SCALAIRE ENTRE LA NORME DE LA SURFACE ET DR
!       SEUIL  : SEUIL DE LA SURFACE DE CHARGE ANTERIEURE
!   -------------------------------------------------------------------
#include "asterfort/hujprj.h"
    integer(kind=8) :: ndt, ndi, i, k
    real(kind=8) :: mater(22, 2), sig(6), vin(*)
    real(kind=8) :: b, pco, beta, pc, epsvpm, ptrac
    real(kind=8) :: un, zero
    real(kind=8) :: p, q, m, phi, degr, sigd(3)
    real(kind=8) :: posf(3), ref(2), norm(2), pst, tole1
!
    parameter(degr=0.0174532925199d0)
!
    common/tdim/ndt, ndi
!
    data un, zero, tole1/1.d0, 0.d0, 1.d-7/
!
    b = mater(4, 2)
    pco = mater(7, 2)
    beta = mater(2, 2)
    epsvpm = vin(23)
    phi = mater(5, 2)
    m = sin(degr*phi)
    ptrac = mater(21, 2)
!
    pc = pco*exp(-beta*epsvpm)
!
    call hujprj(k, sig, sigd, p, q)
!
    p = p-ptrac
!
    do i = 1, 3
        if (q .gt. tole1) then
            posf(i) = sigd(i)/(m*p*(un-b*log(p/pc)))
        else
            posf(i) = zero
        end if
    end do
    norm(1) = vin(4*k+7)
    norm(2) = vin(4*k+8)
    ref(1) = vin(4*k+5)
    ref(2) = vin(4*k+6)
!
    pst = 2.d0*norm(1)*(posf(1)-ref(1))+norm(2)*(posf(3)-ref(2))
!
end subroutine
