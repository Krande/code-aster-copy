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

! --------------------------------------------------------------------------------------------------
!  Damage criterion: compute the damage measure chi and its derivative
!  chi = < eigmax(tns) >
! --------------------------------------------------------------------------------------------------

module endo_crit_tc_module

    use tenseur_dime_module, only: rs

    implicit none
    private
    public:: CRITERION, Init, Derivative

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcvalp.h"
#include "asterfort/lcesme.h"

! --------------------------------------------------------------------

    ! Damage criterion
    type CRITERION
        integer(kind=8)                 :: ndimsi
        real(kind=8)            :: p
        real(kind=8)            :: prec
        real(kind=8), allocatable:: tns(:)
        real(kind=8)            :: eig(3)
        real(kind=8)            :: chi
    end type CRITERION

contains

    function Init(p, tns, prec) result(self)

        implicit none
        real(kind=8), intent(in)    :: p
        real(kind=8), intent(in)    :: tns(:)
        real(kind=8), intent(in)    :: prec
        type(CRITERION)             :: self
! --------------------------------------------------------------------------------------------------
!  OBJECT CREATION AND INITIALISATION
! p         Material characterics of the damage surface -> regularisation parameter
! tns       driving tensor
! prec      Accuracy for derivative computation
! --------------------------------------------------------------------------------------------------

        allocate (self%tns(size(tns)))

        self%p = p
        self%ndimsi = size(tns)
        self%prec = prec
        self%tns = tns
        call ComputeChi(self)

    end function Init

    subroutine ComputeChi(self)
        implicit none
        type(CRITERION)        :: self
! --------------------------------------------------------------------------------------------------
! Compute the damage measure
! --------------------------------------------------------------------------------------------------
        real(kind=8):: poseigmax
! --------------------------------------------------------------------------------------------------

        ! eigenvalues
        call lcvalp(rs(6, self%tns), self%eig)
        poseigmax = maxval(max(self%eig, 0.d0))

        ! Positive part of the maximal eigenvalue (with control to avoid overflow)
        if (poseigmax .le. 1.d0) then
            self%chi = sum(max(self%eig, 0.d0)**self%p)**(1.d0/self%p)
        else
            self%chi = poseigmax*sum(max(self%eig/poseigmax, 0.d0)**self%p)**(1.d0/self%p)
        end if

    end subroutine ComputeChi

    function Derivative(self) result(dtns_chi)
        implicit none
        type(CRITERION)        :: self
        real(kind=8)           :: dtns_chi(1:self%ndimsi)
! --------------------------------------------------------------------------------------------------
! Compute the derivative of chi with respect to the driving tensor
! --------------------------------------------------------------------------------------------------
        real(kind=8), parameter:: safe = 1.d-2
! --------------------------------------------------------------------------------------------------
        real(kind=8):: para(1), fct_6(6), drv_66(6, 6)
! --------------------------------------------------------------------------------------------------

        ! Probleme potentiel dans le cas non local ?
        ASSERT(self%chi .gt. 0.d0)

        ! Tensorial derivative (with control to avoid overflow)
        para(1) = self%p
        call lcesme(rs(6, self%tns/self%chi), self%eig/self%chi, para, power_pm1, &
                    self%prec*safe, fct_6, drv_66)
        dtns_chi = fct_6(1:self%ndimsi)

    end function Derivative

    subroutine power_pm1(x, para, fct, der)
        implicit none

        real(kind=8), intent(in) :: x, para(:)
        real(kind=8), intent(out):: fct, der
! --------------------------------------------------------------------------------------------------
! Function x -> <x>**(p-1) pour le calcul de la fonction tensorielle correspondante
! Attention: derivee fixee a zero pour eviter les problemes de calcul (elle n'est pas utilisee)
! x      argument
! para   additional parameters
!        para(1) = p (exponent)
! fct    f(x)
! der    f'(x)
! --------------------------------------------------------------------------------------------------
        real(kind=8) :: p
! --------------------------------------------------------------------------------------------------

        ASSERT(size(para) .eq. 1)
        p = para(1)

        if (x .gt. 0.d0) then
            fct = x**(p-1)
            der = 0.0
!        der = (p-1)*(x**(p-2))
        else
            fct = 0.d0
            der = 0.d0
        end if

    end subroutine power_pm1

end module endo_crit_tc_module
