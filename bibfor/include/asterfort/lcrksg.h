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
    subroutine lcrksg(rela_comp, nvi, vinf, fd, df,&
                      nmat, coefl, sigi)
        integer(kind=8) :: nmat
        character(len=16) :: rela_comp
        integer(kind=8) :: nvi
        real(kind=8) :: vinf(*)
        real(kind=8) :: fd(9)
        real(kind=8) :: df(9)
        real(kind=8) :: coefl(nmat)
        real(kind=8) :: sigi(6)
    end subroutine lcrksg
end interface
