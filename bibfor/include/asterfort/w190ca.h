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
    subroutine w190ca(modele, carele, chmar1, chefge, chamfer, chefge0, chmar2)
        character(len=8) :: modele
        character(len=8) :: carele
        character(len=19) :: chmar1
        character(len=19) :: chefge
        character(len=19) :: chamfer
        character(len=19) :: chefge0
        character(len=19) :: chmar2
    end subroutine w190ca
end interface
