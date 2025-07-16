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

module endo_rigi_unil_module

    use tenseur_dime_module, only: kron, rs

    implicit none
    private
    public:: MATERIAL, UNILATERAL, Init, ComputeEnergy, ComputeStress, ComputeStress_eig, &
             ComputeStiffness

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcvalp.h"
#include "asterfort/lcesme.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics

    type MATERIAL
        real(kind=8) :: lambda, deuxmu, regbet
    end type MATERIAL

    ! Unilateral treatment
    type UNILATERAL
        type(MATERIAL)                           :: mat
        integer(kind=8)                                  :: ndimsi
        real(kind=8), dimension(3)             :: eps_eig, sigall_eig, sigpos_eig, signeg_eig
        real(kind=8)                             :: treps
        real(kind=8)                             :: unitr
        real(kind=8)                             :: dertr
        real(kind=8), dimension(:), allocatable    :: eps, unieps, sigall, sigpos, signeg
        real(kind=8), dimension(:, :), allocatable  :: dereps
    end type UNILATERAL

contains

! ==================================================================================================
!  OBJECT CREATION AND INITIALISATION
! ==================================================================================================

    function Init(mat, eps, prece) result(self)

        implicit none
        type(MATERIAL), intent(in)  :: mat
        real(kind=8), intent(in)    :: eps(:)
        real(kind=8), intent(in)    :: prece
        type(UNILATERAL)            :: self
! --------------------------------------------------------------------------------------------------
! mat       Material characterics
! eps       Current strain state
! prece     relative accuracy with respect to the strain components
! --------------------------------------------------------------------------------------------------
        integer(kind=8):: i
        real(kind=8):: para(2), rdum
        real(kind=8), dimension(size(eps)):: kr
        real(kind=8), dimension(3)::eps_eig, negeps_eig
        real(kind=8), dimension(6)::unieps_6
        real(kind=8), dimension(6, 6)::dereps_66
! --------------------------------------------------------------------------------------------------
        real(kind=8), parameter::safe = 1.d2
! --------------------------------------------------------------------------------------------------

        ! Size allocation
        self%ndimsi = size(eps)
        allocate (self%eps(self%ndimsi))
        allocate (self%unieps(self%ndimsi))
        allocate (self%dereps(self%ndimsi, self%ndimsi))

        ! Initialisation
        kr = kron(self%ndimsi)
        para(1) = mat%regbet
        para(2) = 0.d0

        self%mat = mat
        self%eps = eps
        self%treps = sum(eps(1:3))

        ! Trace term
        call NegPart(self%treps, para, self%unitr, self%dertr)

        ! Eigenvalues in decreasing order
        call lcvalp(rs(6, eps), eps_eig)
        self%eps_eig = eps_eig

        ! Value and derivative of the tensorial unilateral function (negative part)
        call lcesme(rs(6, eps), eps_eig, para, NegPart, prece/safe, unieps_6, dereps_66)
        self%unieps = unieps_6(1:self%ndimsi)
        self%dereps = dereps_66(1:self%ndimsi, 1:self%ndimsi)

        ! Principal stresses
        do i = 1, 3
            call NegPart(self%eps_eig(i), para, negeps_eig(i), rdum)
        end do
        self%sigall_eig = self%mat%lambda*self%treps+self%mat%deuxmu*self%eps_eig
        self%signeg_eig = self%mat%lambda*self%unitr+self%mat%deuxmu*negeps_eig
        self%sigpos_eig = self%sigall_eig-self%signeg_eig

        ! Stresses
        self%sigall = self%mat%lambda*self%treps*kr+self%mat%deuxmu*self%eps
        self%signeg = self%mat%lambda*self%unitr*kr+self%mat%deuxmu*self%unieps
        self%sigpos = self%sigall-self%signeg

    end function Init

! ==================================================================================================
!  ENERGIE AVEC RESTAURATION DE RIGIDITE
! ==================================================================================================

    subroutine ComputeEnergy(self, wpos, wneg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: wpos, wneg
! --------------------------------------------------------------------------------------------------
! wpos      (regularised) tensile contribution to the energy
! wneg      (regularised) compressive contribution to the energy
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: i
        real(kind=8):: para(1)
        real(kind=8):: trNhs, wall
        real(kind=8), dimension(3):: eigNhs
! --------------------------------------------------------------------------------------------------

!   Initialisation
        para(1) = self%mat%regbet

!   Total energy
        wall = 0.5d0*(self%mat%lambda*self%treps**2+self%mat%deuxmu*sum(self%eps_eig**2))

!   Negative and positive parts
        trNhs = NegHalfSquare(self%treps, para)
        eigNhs = [(NegHalfSquare(self%eps_eig(i), para), i=1, size(self%eps_eig))]
        wneg = self%mat%lambda*trNhs+self%mat%deuxmu*sum(eigNhs)
        wpos = wall-wneg

    end subroutine ComputeEnergy

! ==================================================================================================
!  CONTRAINTES AVEC RESTAURATION DE RIGIDITE
! ==================================================================================================

    subroutine ComputeStress(self, sigpos, signeg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: sigpos(:), signeg(:)
! --------------------------------------------------------------------------------------------------
! sigpos    contrainte positive (traction)
! signeg    contrainte negative (compression)
! --------------------------------------------------------------------------------------------------

        sigpos = self%sigpos
        signeg = self%signeg

    end subroutine ComputeStress

! ==================================================================================================
!  CONTRAINTES AVEC RESTAURATION DE RIGIDITE
! ==================================================================================================

    subroutine ComputeStress_eig(self, sigpos_eig, signeg_eig)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: sigpos_eig(3), signeg_eig(3)
! --------------------------------------------------------------------------------------------------
! sigpos_eig    contrainte principales de traction (ordre decroissant)
! signeg_eig    contrainte principales de compression (ordre devroissant)
! --------------------------------------------------------------------------------------------------

        sigpos_eig = self%sigpos_eig
        signeg_eig = self%signeg_eig

    end subroutine ComputeStress_eig

! ==================================================================================================
!  RIGIDITE UNILATERALES
! ==================================================================================================

    subroutine ComputeStiffness(self, de_spos, de_sneg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: de_spos(:, :), de_sneg(:, :)
! --------------------------------------------------------------------------------------------------
! de_spos   derivee de la contrainte sigpos / deformation
! de_sneg   derivee de la contrainte signeg / deformation
! --------------------------------------------------------------------------------------------------
        integer(kind=8):: i
        real(kind=8), dimension(self%ndimsi, self%ndimsi):: de_sall
        real(kind=8), dimension(self%ndimsi):: kr
! --------------------------------------------------------------------------------------------------

!   Initialisation
        kr = kron(self%ndimsi)

        ! Matrice totale (Hooke)
        de_sall = 0.d0
        de_sall(1:3, 1:3) = self%mat%lambda
        do i = 1, self%ndimsi
            de_sall(i, i) = de_sall(i, i)+self%mat%deuxmu
        end do

        ! Matrices negative et positive
        de_sneg = self%mat%deuxmu*self%dereps
        de_sneg(1:3, 1:3) = de_sneg(1:3, 1:3)+self%mat%lambda*self%dertr
        de_spos = de_sall-de_sneg

    end subroutine ComputeStiffness

! ==================================================================================================
!   SMOOTHED NEGATIVE HALF SQUARE FUNCTION
!   f(x) approximate 0.5 * <-x>**2
!     f(x) = 0.5 * x**2 * exp(beta/x)) if x<0
!     f(x) = 0                             if x>0
! ==================================================================================================

    function NegHalfSquare(x, p) result(nhs)
        implicit none

        real(kind=8), intent(in):: x, p(:)
        real(kind=8)           :: nhs
! --------------------------------------------------------------------------------------------------
! x:   array argument (the function is applied to each x(i))
! p:   additional parameters
!        p(1) = beta (smoothing parameter)
! --------------------------------------------------------------------------------------------------
        real(kind=8) :: beta
! --------------------------------------------------------------------------------------------------

        ASSERT(size(p) .eq. 1)
        beta = p(1)
        if (x .ge. -1.d-3*beta) then
            nhs = 0.d0
        else
            nhs = 0.5d0*x**2*exp(beta/x)
        end if

    end function NegHalfSquare

! ==================================================================================================
!   SMOOTHED UNILATERAL FUNCTION AND ITS DERIVATIVE
!     f(x) = (x - 0.5*beta) * exp(beta/x)) if x<0
!     f(x) = 0                                  if x>0
! ==================================================================================================

    subroutine NegPart(x, p, fct, der)
        implicit none

        real(kind=8), intent(in) :: x, p(:)
        real(kind=8), intent(out):: fct, der
! --------------------------------------------------------------------------------------------------
! x:   argument
! p:   additional parameters
!        p(1) = beta (smoothing parameter)
!        p(2) = 0 si secant, 1 si tangent
! fct: f(x)
! der: f'(x) ou f(x)/x si secant
! --------------------------------------------------------------------------------------------------
        aster_logical:: elas
        real(kind=8) :: beta, u
! --------------------------------------------------------------------------------------------------

        ASSERT(size(p) .eq. 2)
        beta = p(1)
        elas = p(2) .gt. 0.5d0

        if (x .ge. -1.d-3*beta) then
            fct = 0
            der = 0
        else
            u = -beta/x
            fct = (x-0.5d0*beta)*exp(beta/x)
            der = merge(fct/x, (1+u*(1+0.5d0*u))*exp(-u), elas)
        end if

    end subroutine NegPart

end module endo_rigi_unil_module
