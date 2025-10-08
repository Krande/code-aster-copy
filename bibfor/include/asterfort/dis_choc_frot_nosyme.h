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
interface
    subroutine dis_choc_frot_nosyme(for_discret, icodma, ulp, xg, klv, &
                             dpe, varmo, force, varpl)
        use te0047_type
        type(te0047_dscr), intent(in) :: for_discret
        integer(kind=8)         :: icodma
        real(kind=8)    :: ulp(*), klv(*), xg(*), varmo(*), varpl(*), force(*), dpe(*)
    end subroutine dis_choc_frot_nosyme
end interface
