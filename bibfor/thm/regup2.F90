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

subroutine regup2(x0, y0, y0p, y1, &
                  a, b, c)
!
    implicit none
!
    real(kind=8), intent(in) :: x0, y0, y0p, y1
    real(kind=8), intent(out) :: a, b, c
!
! --------------------------------------------------------------------------------------------------
!
! THM - Permeability
!
! Regularization with quadratic polynom - Get coefficients of polynom
!
! --------------------------------------------------------------------------------------------------
!
! For polynom p=ax²+bx+c with p(x0) = y0, dp_dx(x0) = y0p and p(1) = y1
!
! --------------------------------------------------------------------------------------------------
!
    a = (y1-y0+y0p*(x0-1.d0))/((1.d0-x0)**2)
    b = -2.d0*a*x0+y0p
    c = y1-a-b
!
end subroutine
