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
    subroutine wpermo(lmasse, lraide, lamor, nbprop, vecp,&
                      fr, am, excl, omecor, ernorm)
        integer(kind=8) :: lmasse
        integer(kind=8) :: lraide
        integer(kind=8) :: lamor
        integer(kind=8) :: nbprop
        complex(kind=8) :: vecp(*)
        real(kind=8) :: fr(*)
        real(kind=8) :: am(*)
        integer(kind=8) :: excl(*)
        real(kind=8) :: omecor
        real(kind=8) :: ernorm(*)
    end subroutine wpermo
end interface
