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

subroutine lcesus(x, p, val, der)
    implicit none

    real(kind=8), intent(in) :: x, p(:)
    real(kind=8), intent(out):: val, der
! --------------------------------------------------------------------------------------------------
!   smoothed unilateral function and its pseudo (secant) derivative:
!     f(x) = (x - 0.5/gamma) * exp(1/(gamma*x)) if x<0
!     f(x) = 0                                  if x>0
! --------------------------------------------------------------------------------------------------
! x:   argument
! p:   dummy argument
! val: f(x)
! der: f'(x)
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
    common/lcee/lambda, deuxmu, troisk, gamma, rigmin, pc, pr, epsth
! --------------------------------------------------------------------------------------------------
    if (x*gamma .ge. -1.d-3) then
        val = 0
        der = 0
    else
        val = (x-0.5d0/gamma)*exp(1/(gamma*x))
        der = val/x
    end if

end subroutine lcesus
