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
subroutine separ_RI_elas_cp(elas_id, nu, g, nui, gi, &
                            e1, e2, &
                            nu12, &
                            e1i, e2i, &
                            nu12i, &
                            hr, hi)
!
    implicit none
!
    integer(kind=8), intent(in) :: elas_id
    real(kind=8), intent(in) :: nu, g, e1, e2
    real(kind=8), intent(in) :: nu12
    real(kind=8), intent(in) :: nui, gi, e1i, e2i
    real(kind=8), intent(in) :: nu12i
    real(kind=8), intent(out) :: hr(3), hi(3)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Separate elastic modulus for the computation of imaginary and real part of the rigidity
!
! --------------------------------------------------------------------------------------------------
!
! In  elas_id          : Type of elasticity or viscoelasticity
!                 1 - Isotropic elasticity
!                 2 - Orthotropic elasticity
!                 3 - Transverse isotropic elasticity
!                 4 - Isotropic elasticity
!                 5 - Orthotropic elasticity
!                 6 - Transverse isotropic elasticity
! In nu                : real Poisson ratio (isotropic)
! In nui               : imaginary Poisson ratio (isotropic)
! In g                 : real shear ratio (isotropic/Transverse isotropic)
! In gi                : imaginary shear ratio (isotropic)
! IN e1                : real Young modulus - Direction 1 (Orthotropic/Transverse isotropic)
! In e2                : real Young modulus - Direction 2 (Orthotropic)
! In e1i               : imaginary Young modulus - Direction 1 (Orthotropic/Transverse isotropic)
! In e2i               : imaginary Young modulus - Direction 2 (Orthotropic)
! In nu12              : real Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! In nu12i             : imaginary Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! hr                   : real part of Hook matrix component
! hi                   : imaginary part of Hook matrix component
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nu21, delta, c1
    complex(kind=8) :: e1c, e2c, nu12c, nu21c
    complex(kind=8) :: nuc, Gc, deltac, c1c
    real(kind=8), parameter :: undemi = 0.5d0
    real(kind=8), parameter :: un = 1.d0
!
! --------------------------------------------------------------------------------------------------
!
    if (elas_id .eq. 1) then
!
        hr(1) = 2.d0*g/(1.d0-nu)
        hr(2) = nu*hr(1)

!
    elseif (elas_id .eq. 4) then
!
        nuc = dcmplx(nu, nui)
        Gc = dcmplx(g, gi)
!
        hr(1) = real(2.d0*Gc/(1.d0-nuc))
        hr(2) = real(2.d0*Gc*nuc/(1.d0-nuc))
!
        hi(1) = aimag(2.d0*Gc/(1.d0-nuc))
        hi(2) = aimag(2.d0*Gc*nuc/(1.d0-nuc))
!
    elseif (elas_id .eq. 2) then
!
        nu21 = e2*nu12/e1
        delta = un-nu12*nu21
        hr(1) = e1/delta
        hr(2) = nu12*e2/delta
        hr(3) = e2/delta
!
    elseif (elas_id .eq. 5) then
!
        e1c = dcmplx(e1, e1i)
        e2c = dcmplx(e2, e2i)
        nu12c = dcmplx(nu12, nu12i)
!
        nu21c = e2c*nu12c/e1c
        deltac = un-nu12c*nu21c
        hr(1) = real(e1c/deltac)
        hr(2) = real(nu12c*e2c/deltac)
        hr(3) = real(e2c/deltac)
        hi(1) = aimag(e1c/deltac)
        hi(2) = aimag(nu12c*e2c/deltac)
        hi(3) = aimag(e2c/deltac)
!
    elseif (elas_id .eq. 3) then
!
        c1 = e1/(un+nu12)
        delta = un-nu12*nu12
        hr(1) = e1/delta
        hr(2) = nu12*hr(1)
        hr(3) = undemi*c1
!
    elseif (elas_id .eq. 6) then
!
        e1c = dcmplx(e1, e1i)
        nu12c = dcmplx(nu12, nu12i)
!
        c1c = e1c/(un+nu12c)
        deltac = un-nu12c*nu12c
        hr(1) = real(e1c/deltac)
        hr(2) = real(nu12*e1c/deltac)
        hr(3) = real(undemi*c1c)
        hi(1) = aimag(e1c/deltac)
        hi(2) = aimag(nu12*e1c/deltac)
        hi(3) = aimag(undemi*c1c)
!
    end if
!
end subroutine
