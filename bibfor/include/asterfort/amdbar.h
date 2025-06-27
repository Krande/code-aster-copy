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
    subroutine amdbar(n, pe, iw, len, iwlen,&
                      pfree, nv, next, last, head,&
                      elen, degree, ncmpa, w, iovflo)
        integer(kind=8) :: iwlen
        integer(kind=8) :: n
        integer(kind=8) :: pe(n)
        integer(kind=8) :: iw(iwlen)
        integer(kind=8) :: len(n)
        integer(kind=8) :: pfree
        integer(kind=8) :: nv(n)
        integer(kind=8) :: next(n)
        integer(kind=8) :: last(n)
        integer(kind=8) :: head(n)
        integer(kind=8) :: elen(n)
        integer(kind=8) :: degree(n)
        integer(kind=8) :: ncmpa
        integer(kind=8) :: w(n)
        integer(kind=8) :: iovflo
    end subroutine amdbar
end interface
