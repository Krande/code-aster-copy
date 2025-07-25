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
    subroutine ptforp(itype, option, nomte, a, a2, &
                      xl, ist, nno, ncf, pgl, fer,  fei)
        integer(kind=8)             :: itype
        character(len=*)    :: option
        character(len=*)    :: nomte
        real(kind=8)        :: a
        real(kind=8)        :: a2
        real(kind=8)        :: xl
        integer(kind=8)             :: ist
        integer(kind=8)             :: nno
        integer(kind=8)             :: ncf
        real(kind=8)        :: pgl(3, 3)
        real(kind=8)        :: fer(*)
        real(kind=8)        :: fei(*)
    end subroutine ptforp
end interface
