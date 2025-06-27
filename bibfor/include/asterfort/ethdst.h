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
    subroutine ethdst(fami, nno, ndim, nbsig, npg,&
                      ipoids, ivf, idfde, xyz, depl,&
                      instan, angl_naut, mater, option, enthth)
        character(len=*) :: fami
        integer(kind=8) :: nno
        integer(kind=8) :: ndim
        integer(kind=8) :: nbsig
        integer(kind=8) :: npg
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        real(kind=8) :: xyz(*)
        real(kind=8) :: depl(*)
        real(kind=8) :: instan
        real(kind=8) :: angl_naut(3)
        integer(kind=8) :: mater
        character(len=16) :: option
        real(kind=8) :: enthth
    end subroutine ethdst
end interface
