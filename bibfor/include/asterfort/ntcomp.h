! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
#include "asterf_types.h"
!
interface
    subroutine ntcomp(rela_name, ndim, temp, dtemp, coorpg, aniso, ifon, fluxglo, Kglo, dfluxglo)
        character(len=16), intent(in) :: rela_name
        integer(kind=8), intent(in) :: ndim, ifon(6)
        real(kind=8), intent(in) :: temp, dtemp(3), coorpg(3)
        aster_logical, intent(in) :: aniso
        real(kind=8), intent(out) :: fluxglo(3)
        real(kind=8), intent(out) :: Kglo(3, 3), dfluxglo(3)
    end subroutine ntcomp
end interface
