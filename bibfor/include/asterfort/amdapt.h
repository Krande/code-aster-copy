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
    subroutine amdapt(neq, nbnd, nbsn, pe, nv,&
                      invp, parent, supnd, adress, lgind,&
                      fctnzs, fctops, llist, nnv)
        integer(kind=8) :: neq
        integer(kind=8) :: nbnd
        integer(kind=8) :: nbsn
        integer(kind=8) :: pe(neq+1)
        integer(kind=8) :: nv(neq)
        integer(kind=8) :: invp(neq)
        integer(kind=8) :: parent(*)
        integer(kind=8) :: supnd(neq)
        integer(kind=8) :: adress(*)
        integer(kind=8) :: lgind
        integer(kind=8) :: fctnzs
        real(kind=8) :: fctops
        integer(kind=8) :: llist(neq)
        integer(kind=8) :: nnv(neq)
    end subroutine amdapt
end interface
