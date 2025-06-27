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
    subroutine crvloc(dim, adcom0, iatyma, jconnex0, jconnexc, vgeloc,&
                      nvtot, nvoima, nscoma, touvoi)
        integer(kind=8) :: nscoma
        integer(kind=8) :: nvoima
        integer(kind=8), intent(in) :: dim
        integer(kind=8) :: adcom0
        integer(kind=8) :: iatyma
        integer(kind=8), intent(in) :: jconnex0
        integer(kind=8), intent(in) :: jconnexc
        integer(kind=8) :: vgeloc(*)
        integer(kind=8) :: nvtot
        integer(kind=8) :: touvoi(1:nvoima, 1:nscoma+2)
    end subroutine crvloc
end interface
