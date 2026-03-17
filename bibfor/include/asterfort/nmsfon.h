! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
interface
    subroutine nmsfon(refe, ndim, nno, npg, nddl, &
                      geomi, vff, idff, iw,  sief, fint)

        integer(kind=8), intent(in)          :: ndim, nno, npg, nddl
        integer(kind=8), intent(in)          :: iw, idff
        real(kind=8), intent(in)     :: geomi(ndim, nno), vff(nno, npg)
        real(kind=8), intent(in)     :: sief(4*ndim, npg)
        aster_logical, intent(in)    :: refe
        real(kind=8), intent(out)    :: fint(nddl)
    end subroutine nmsfon
end interface
