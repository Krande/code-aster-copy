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
!
interface
    subroutine pmppr(amat, na1, na2, ka, bmat,&
                     nb1, nb2, kb, cmat, nc1,&
                     nc2)
        integer(kind=8), intent(in) :: ka
        integer(kind=8), intent(in) :: kb
        integer(kind=8), intent(in) :: na1
        integer(kind=8), intent(in) :: na2
        integer(kind=8), intent(in) :: nb1
        integer(kind=8), intent(in) :: nb2
        integer(kind=8), intent(in) :: nc1
        integer(kind=8), intent(in) :: nc2
        real(kind=8), intent(in) :: amat(na1, na2)
        real(kind=8), intent(in) :: bmat(nb1, nb2)
        real(kind=8), intent(out) :: cmat(nc1, nc2)
    end subroutine pmppr
end interface
