! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
    public:: MATERIAL, UNILATERAL, Init, ComputeEnergy, ComputeStress, ComputeStiffness

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcvalp.h"
#include "asterfort/lcesme.h"
#include "asterc/r8vide.h"

! --------------------------------------------------------------------

    ! Material characteristics

    type MATERIAL
        real(kind=8) :: lambda, deuxmu, regam
    end type MATERIAL

    ! Unilateral treatment
    type UNILATERAL
        type(MATERIAL)                           :: mat
        integer                                  :: ndimsi
        real(kind=8)                             :: eigeps(3)
        real(kind=8)                             :: treps
        real(kind=8)                             :: unitr
        real(kind=8)                             :: dertr
        real(kind=8), dimension(:), allocatable    :: eps
        real(kind=8), dimension(:), allocatable    :: unieps
        real(kind=8), dimension(:, :), allocatable  :: dereps
    end type UNILATERAL

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(mat, eps, prece) result(self)

        implicit none
        type(MATERIAL), intent(in)   :: mat
        real(kind=8), intent(in)    :: eps(:)
        real(kind=8), intent(in)    :: prece
        type(UNILATERAL)            :: self
! --------------------------------------------------------------------------------------------------
! mat       Material characterics
! eps       Current strain state
! prece     relative accuracy with respect to the strain components
! --------------------------------------------------------------------------------------------------
        real(kind=8):: para(2)
        real(kind=8), dimension(3)::eigeps
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
        para(1) = mat%regam
        para(2) = 0.d0

        self%mat = mat
        self%eps = eps
        self%treps = sum(eps(1:3))

        ! Trace term
        call NegPart(self%treps, para, self%unitr, self%dertr)

        ! Eigenvalues
        call lcvalp(rs(6, eps), eigeps)
        self%eigeps = eigeps

        ! Value and derivative of the tensorial unilateral function (negative part)
        call lcesme(rs(6, eps), eigeps, para, NegPart, prece/safe, unieps_6, dereps_66)
        self%unieps = unieps_6(1:self%ndimsi)
        self%dereps = dereps_66(1:self%ndimsi, 1:self%ndimsi)

    end function Init

! =====================================================================
!  ENERGIE AVEC RESTAURATION DE RIGIDITE
! =====================================================================

    subroutine ComputeEnergy(self, wpos, wneg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: wpos, wneg
! ---------------------------------------------------------------------
! wpos      (regularised) tensile contribution to the energy
! wneg      (regularised) compressive contribution to the energy
! ---------------------------------------------------------------------
        integer     :: i
        real(kind=8):: para(1)
        real(kind=8):: trNhs, wall
        real(kind=8), dimension(3):: eigNhs
! ---------------------------------------------------------------------

!   Initialisation
        para(1) = self%mat%regam

!   Total energy
        wall = 0.5d0*(self%mat%lambda*self%treps**2+self%mat%deuxmu*sum(self%eigeps**2))

!   Negative and positive parts
        trNhs = NegHalfSquare(self%treps, para)
        eigNhs = [(NegHalfSquare(self%eigeps(i), para), i=1, size(self%eigeps))]
        wneg = self%mat%lambda*trNhs+self%mat%deuxmu*sum(eigNhs)
        wpos = wall-wneg

    end subroutine ComputeEnergy

! =====================================================================
!  CONTRAINTES ET RIGIDITE AVEC RESTAURATION DE RIGIDITE
! =====================================================================

    subroutine ComputeStress(self, sigpos, signeg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: sigpos(:), signeg(:)
! ---------------------------------------------------------------------
! sigpos    contrainte positive (traction)
! signeg    contrainte negative (compression)
! ---------------------------------------------------------------------
        real(kind=8), dimension(self%ndimsi):: sigall, kr
! ---------------------------------------------------------------------

!   Initialisation
        kr = kron(self%ndimsi)

!   Stresses
        sigall = self%mat%lambda*self%treps*kr+self%mat%deuxmu*self%eps
        signeg = self%mat%lambda*self%unitr*kr+self%mat%deuxmu*self%unieps
        sigpos = sigall-signeg

    end subroutine ComputeStress

! =====================================================================
!  RIGIDITE UNILATERALES
! =====================================================================

    subroutine ComputeStiffness(self, de_spos, de_sneg)
        implicit none

        type(UNILATERAL), intent(in):: self
        real(kind=8), intent(out)   :: de_spos(:, :), de_sneg(:, :)
! ---------------------------------------------------------------------
! de_spos   derivee de la contrainte sigpos / deformation
! de_sneg   derivee de la contrainte signeg / deformation
! ---------------------------------------------------------------------
        integer:: i
        real(kind=8), dimension(self%ndimsi, self%ndimsi):: de_sall
        real(kind=8), dimension(self%ndimsi):: kr
! ---------------------------------------------------------------------

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

! =====================================================================
!   SMOOTHED NEGATIVE HALF SQUARE FUNCTION
!   f(x) approximate 0.5 * <-x>**2
!     f(x) = 0.5 * x**2 * exp(1/(gamma*x)) if x<0
!     f(x) = 0                             if x>0
! =====================================================================

    function NegHalfSquare(x, p) result(nhs)
        implicit none

        real(kind=8), intent(in):: x, p(:)
        real(kind=8)           :: nhs
! ---------------------------------------------------------------------
! x:   array argument (the function is applied to each x(i))
! p:   additional parameters
!        p(1) = gamma (smoothing parameter)
! ---------------------------------------------------------------------
        real(kind=8) :: ga
! ---------------------------------------------------------------------

        ASSERT(size(p) .eq. 1)
        ga = p(1)
        if (ga*x .ge. -1.d-3) then
            nhs = 0.d0
        else
            nhs = 0.5d0*x**2*exp(1/(ga*x))
        end if

    end function NegHalfSquare

! =====================================================================
!   SMOOTHED UNILATERAL FUNCTION AND ITS DERIVATIVE
!     f(x) = (x - 0.5/gamma) * exp(1/(gamma*x)) if x<0
!     f(x) = 0                                  if x>0
! =====================================================================

    subroutine NegPart(x, p, fct, der)
        implicit none

        real(kind=8), intent(in) :: x, p(:)
        real(kind=8), intent(out):: fct, der
! ---------------------------------------------------------------------
! x:   argument
! p:   additional parameters
!        p(1) = gamma (smoothing parameter)
!        p(2) = 0 si secant, 1 si tangent
! fct: f(x)
! der: f'(x) ou f(x)/x si secant
! ---------------------------------------------------------------------
        aster_logical:: elas
        real(kind=8) :: regam
! ---------------------------------------------------------------------

        ASSERT(size(p) .eq. 2)
        regam = p(1)
        elas = p(2) .gt. 0.5d0

        if (x*regam .ge. -1.d-3) then
            fct = 0
            der = 0
        else
            fct = (x-0.5d0/regam)*exp(1/(regam*x))
            der = merge(fct/x, (1-(regam*x-0.5d0)/(regam*x)**2)*exp(1/(regam*x)), elas)
        end if

    end subroutine NegPart

end module endo_rigi_unil_module
