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
    subroutine xvechp(ndim, elrefp, nnop, igeom, itemp,&
                      itps, ihechp, jptint, jcface,&
                      jlonch, jlst, jbasec, nfh, nfe,&
                      fonree, ivectt, heavn)
        integer(kind=8) :: nfe
        integer(kind=8) :: nfh
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        integer(kind=8) :: igeom
        integer(kind=8) :: itemp
        integer(kind=8) :: itps
        integer(kind=8) :: ihechp
        integer(kind=8) :: jptint
        integer(kind=8) :: jcface
        integer(kind=8) :: jlonch
        integer(kind=8) :: jlst
        integer(kind=8) :: jbasec
        integer(kind=8) :: heavn(27,5)
        character(len=4) :: fonree
        integer(kind=8) :: ivectt
    end subroutine xvechp
end interface
