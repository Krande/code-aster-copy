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
    subroutine arltem(ndim  ,nomte, &
                      nns   ,jcoors, &
                      npgs  ,ivfs  ,idfdes,ipoids, &
                      elref1,ndml1   ,jcoor1, &
                      elref2,ndml2   ,jcoor2, &
                      mcpln1,mcpln2)
        integer(kind=8) :: ndim
        integer(kind=8) :: nns
        integer(kind=8) :: npgs
        integer(kind=8) :: jcoor1
        integer(kind=8) :: jcoor2
        integer(kind=8) :: jcoors
        integer(kind=8) :: ivfs
        integer(kind=8) :: ipoids
        integer(kind=8) :: idfdes
        integer(kind=8) :: ndml1
        integer(kind=8) :: ndml2
        character(len=8) :: elref1
        character(len=8) :: elref2
        character(len=16) :: nomte
        real(kind=8) :: mcpln1(2*ndim*ndml2,ndim*ndml1)
        real(kind=8) :: mcpln2(2*ndim*ndml2,2*ndim*ndml2)
    end subroutine arltem
end interface
