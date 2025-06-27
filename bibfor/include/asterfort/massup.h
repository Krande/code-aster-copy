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
    subroutine massup(option, ndim, dlns, nno, nnos,&
                      mate, phenom, npg, ipoids, idfde,&
                      geom, vff1, imatuu, icodre, igeom,&
                      ivf)
        integer(kind=8) :: npg
        integer(kind=8) :: nno
        integer(kind=8) :: ndim
        character(len=16) :: option
        integer(kind=8) :: dlns
        integer(kind=8) :: nnos
        integer(kind=8) :: mate
        character(len=16) :: phenom
        integer(kind=8) :: ipoids
        integer(kind=8) :: idfde
        real(kind=8) :: geom(ndim, nno)
        real(kind=8) :: vff1(nno, npg)
        integer(kind=8) :: imatuu
        integer(kind=8) :: icodre(1)
        integer(kind=8) :: igeom
        integer(kind=8) :: ivf
    end subroutine massup
end interface
