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
subroutine separ_RI_elas_3D(elas_id, nu, g, nui, gi, &
                            e1, e2, e3, &
                            nu12, nu13, nu23, &
                            e1i, e2i, e3i, &
                            nu12i, nu13i, nu23i, &
                            hr, hi)
!
    implicit none
!
    integer(kind=8), intent(in) :: elas_id
    real(kind=8), intent(in) :: nu, g, e1, e2, e3
    real(kind=8), intent(in) :: nu12, nu13, nu23
    real(kind=8), intent(in) :: nui, gi, e1i, e2i, e3i
    real(kind=8), intent(in) :: nu12i, nu13i, nu23i
    real(kind=8), intent(out) :: hr(6), hi(6)
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
! In e3                : real Young modulus - Direction 3 (Orthotropic/Transverse isotropic)
! In e1i               : imaginary Young modulus - Direction 1 (Orthotropic/Transverse isotropic)
! In e2i               : imaginary Young modulus - Direction 2 (Orthotropic)
! In e3i               : imaginary Young modulus - Direction 3 (Orthotropic/Transverse isotropic)
! In nu12              : real Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! In nu13              : real Poisson ratio - Coupling 1/3 (Orthotropic/Transverse isotropic)
! In nu23              : real Poisson ratio - Coupling 2/3 (Orthotropic)
! In nu12i             : imaginary Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! In nu13i             : imaginary Poisson ratio - Coupling 1/3 (Orthotropic/Transverse isotropic)
! In nu23i             : imaginary Poisson ratio - Coupling 2/3 (Orthotropic)
! hr                   : real part of Hook matrix component
! hi                   : imaginary part of Hook matrix component
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nu21, nu31, nu32, delta, c1
    complex(kind=8) :: e1c, e2c, e3c, nu12c, nu13c, nu23c, nu21c, nu31c, nu32c
    complex(kind=8) :: nuc, Gc, deltac, c1c
    real(kind=8), parameter :: undemi = 0.5d0
    real(kind=8), parameter :: un = 1.d0
    real(kind=8), parameter :: deux = 2.d0
!
! --------------------------------------------------------------------------------------------------
!
    if (elas_id .eq. 1) then
!
        hr(1) = 2.d0*g*(1.d0-nu)/(1.d0-2.d0*nu)
        hr(2) = 2.d0*g*nu/(1.d0-2.d0*nu)
!
    elseif (elas_id .eq. 4) then
!
        nuc = dcmplx(nu, nui)
        Gc = dcmplx(g, gi)
!
        hr(1) = real(2.d0*Gc*(1.d0-nuc)/(1.d0-2.d0*nuc))
        hr(2) = real(2.d0*Gc*nuc/(1.d0-2.d0*nuc))
!
        hi(1) = aimag(2.d0*Gc*(1.d0-nuc)/(1.d0-2.d0*nuc))
        hi(2) = aimag(2.d0*Gc*nuc/(1.d0-2.d0*nuc))
!
    elseif (elas_id .eq. 2) then
!
        nu21 = e2*nu12/e1
        nu31 = e3*nu13/e1
        nu32 = e3*nu23/e2
        delta = un-nu23*nu32-nu31*nu13-nu21*nu12-deux*nu23*nu31*nu12
        hr(1) = (un-nu23*nu32)*e1/delta
        hr(2) = (nu21+nu31*nu23)*e1/delta
        hr(3) = (nu31+nu21*nu32)*e1/delta
        hr(4) = (un-nu13*nu31)*e2/delta
        hr(5) = (nu32+nu31*nu12)*e2/delta
        hr(6) = (un-nu21*nu12)*e3/delta
!
    elseif (elas_id .eq. 5) then
!
        e1c = dcmplx(e1, e1i)
        e2c = dcmplx(e2, e2i)
        e3c = dcmplx(e3, e3i)
        nu12c = dcmplx(nu12, nu12i)
        nu13c = dcmplx(nu13, nu13i)
        nu23c = dcmplx(nu23, nu23i)
!
        nu21c = e2c*nu12c/e1c
        nu31c = e3c*nu13c/e1c
        nu32c = e3c*nu23c/e2c
        deltac = un-nu23c*nu32c-nu31c*nu13c-nu21c*nu12c-deux*nu23c*nu31c*nu12c
        hr(1) = real((un-nu23c*nu32c)*e1c/deltac)
        hr(2) = real((nu21c+nu31c*nu23c)*e1c/deltac)
        hr(3) = real((nu31c+nu21c*nu32c)*e1c/deltac)
        hr(4) = real((un-nu13c*nu31c)*e2c/deltac)
        hr(5) = real((nu32c+nu31c*nu12c)*e2c/deltac)
        hr(6) = real((un-nu21c*nu12c)*e3c/deltac)
        hi(1) = aimag((un-nu23c*nu32c)*e1c/deltac)
        hi(2) = aimag((nu21c+nu31c*nu23c)*e1c/deltac)
        hi(3) = aimag((nu31c+nu21c*nu32c)*e1c/deltac)
        hi(4) = aimag((un-nu13c*nu31c)*e2c/deltac)
        hi(5) = aimag((nu32c+nu31c*nu12c)*e2c/deltac)
        hi(6) = aimag((un-nu21c*nu12c)*e3c/deltac)
!
    elseif (elas_id .eq. 3) then
!
        nu31 = nu13*e3/e1
        c1 = e1/(un+nu12)
        delta = un-nu12-deux*nu13*nu31
        hr(1) = c1*(un-nu13*nu31)/delta
        hr(2) = c1*((un-nu13*nu31)/delta-un)
        hr(3) = e3*nu13/delta
        hr(4) = e3*(un-nu12)/delta
        hr(5) = undemi*c1
!
    elseif (elas_id .eq. 6) then
!
        e1c = dcmplx(e1, e1i)
        e2c = dcmplx(e2, e2i)
        e3c = dcmplx(e3, e3i)
        nu12c = dcmplx(nu12, nu12i)
        nu23c = dcmplx(nu23, nu23i)
        nu13c = dcmplx(nu13, nu13i)
!
        nu31c = nu13c*e3c/e1c
        c1c = e1c/(un+nu12c)
        deltac = un-nu12c-deux*nu13c*nu31c
        hr(1) = real(c1c*(un-nu13c*nu31c)/deltac)
        hr(2) = real(c1c*((un-nu13c*nu31c)/deltac-un))
        hr(3) = real(e3c*nu13c/deltac)
        hr(4) = real(e3c*(un-nu12c)/deltac)
        hr(5) = real(undemi*c1c)
        hi(1) = aimag(c1c*(un-nu13c*nu31c)/deltac)
        hi(2) = aimag(c1c*((un-nu13c*nu31c)/deltac-un))
        hi(3) = aimag(e3c*nu13c/deltac)
        hi(4) = aimag(e3c*(un-nu12c)/deltac)
        hi(5) = aimag(undemi*c1c)
!
    end if
!
end subroutine
