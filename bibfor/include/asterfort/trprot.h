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
    subroutine trprot(model, bamo, tgeom, imodg, iadx,&
                      iady, iadz, isst, iadrp, norm1,&
                      norm2, ndble, num, nu, ma,&
                      mate, mateco, moint, ilires, k, icor)
        character(len=2) :: model
        character(len=8) :: bamo
        real(kind=8) :: tgeom(6)
        integer(kind=8) :: imodg
        integer(kind=8) :: iadx
        integer(kind=8) :: iady
        integer(kind=8) :: iadz
        integer(kind=8) :: isst
        integer(kind=8) :: iadrp
        real(kind=8) :: norm1
        real(kind=8) :: norm2
        integer(kind=8) :: ndble
        character(len=14) :: num
        character(len=14) :: nu
        character(len=8) :: ma
        character(len=*) :: mate, mateco
        character(len=8) :: moint
        integer(kind=8) :: ilires
        integer(kind=8) :: k
        integer(kind=8) :: icor(2)
    end subroutine trprot
end interface
