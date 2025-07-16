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
    subroutine thetapdg(ndim, nno, discr, ff, dfdi, ideg, ilag, ithet, dtdm)
        integer(kind=8), intent(in)       :: ndim
        integer(kind=8), intent(in)       :: nno
        character(len=8)          :: discr
        real(kind=8), intent(in)  :: ff(nno)
        real(kind=8), intent(in)  :: dfdi(nno, ndim)
        integer(kind=8), intent(in)       :: ideg
        integer(kind=8), intent(in)       :: ilag
        integer(kind=8), intent(in)       :: ithet
        real(kind=8), intent(out) :: dtdm(3, 4)
    end subroutine thetapdg
end interface
