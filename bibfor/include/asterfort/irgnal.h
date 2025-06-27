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
    subroutine irgnal(ifi, nbordr, coord, connex, point,&
                      nocmp, nbcmp, numel, nobj, nbel,&
                      cnsc, cnsl, cnsv, partie, jtype,&
                      cnsd)
        integer(kind=8) :: nbcmp
        integer(kind=8) :: ifi
        integer(kind=8) :: nbordr
        real(kind=8) :: coord(*)
        integer(kind=8) :: connex(*)
        integer(kind=8) :: point(*)
        character(len=8) :: nocmp(nbcmp)
        integer(kind=8) :: numel
        character(len=*) :: nobj
        integer(kind=8) :: nbel
        integer(kind=8) :: cnsc(*)
        integer(kind=8) :: cnsl(*)
        integer(kind=8) :: cnsv(*)
        character(len=*) :: partie
        integer(kind=8) :: jtype
        integer(kind=8) :: cnsd(*)
    end subroutine irgnal
end interface
