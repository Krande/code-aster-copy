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
    subroutine xtyhea(nfiss, ifiss, ima, nno, jconx1,&
                      jconx2, jstnl, jstnv, nbheav)
        integer(kind=8) :: nno
        integer(kind=8) :: nfiss
        integer(kind=8) :: ifiss
        integer(kind=8) :: ima
        integer(kind=8) :: jconx1
        integer(kind=8) :: jconx2
        integer(kind=8) :: jstnl(nfiss)
        integer(kind=8) :: jstnv(nfiss)
        integer(kind=8) :: nbheav
    end subroutine xtyhea
end interface
