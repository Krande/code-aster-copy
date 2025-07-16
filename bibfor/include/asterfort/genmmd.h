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
    subroutine genmmd(neqns, neqp1, nadj, xadj, adjncy,&
                      maxint, delta, invp, perm, nbsn,&
                      supnd, adress, parent, gssubs, fctnzs,&
                      fctops, dhead, qsize, llist, marker)
        integer(kind=8) :: nadj
        integer(kind=8) :: neqp1
        integer(kind=8) :: neqns
        integer(kind=8) :: xadj(neqp1)
        integer(kind=8) :: adjncy(nadj)
        integer(kind=8) :: maxint
        integer(kind=8) :: delta
        integer(kind=8) :: invp(neqns)
        integer(kind=8) :: perm(neqns)
        integer(kind=8) :: nbsn
        integer(kind=8) :: supnd(neqp1)
        integer(kind=8) :: adress(neqp1)
        integer(kind=8) :: parent(neqns)
        integer(kind=8) :: gssubs
        integer(kind=8) :: fctnzs
        real(kind=8) :: fctops
        integer(kind=8) :: dhead(neqns)
        integer(kind=8) :: qsize(neqns)
        integer(kind=8) :: llist(neqns)
        integer(kind=8) :: marker(neqns)
    end subroutine genmmd
end interface
