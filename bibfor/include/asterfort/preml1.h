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
#include "asterf_types.h"
!
interface
    subroutine preml1(neq, n2, diag, delg, col,&
                      xadj, adjncy, parent, adress, supnd,&
                      nnz, qsize, llist, suiv, p,&
                      q, invp, perm, lgind, ddlmoy,&
                      nbsn, optnum, lgadjn, nrl, deb,&
                      vois, suit, ier, nec, prno,&
                      deeq, noeud, ddl, invpnd, permnd,&
                      spndnd, xadjd, matgen)
        integer(kind=8) :: lgadjn
        integer(kind=8) :: n2
        integer(kind=8) :: neq
        integer(kind=8) :: diag(0:neq)
        integer(kind=8) :: delg(neq)
        integer(kind=8) :: col(*)
        integer(kind=8) :: xadj(neq+1)
        integer(kind=8) :: adjncy(lgadjn)
        integer(kind=8) :: parent(neq)
        integer(kind=8) :: adress(neq)
        integer(kind=8) :: supnd(neq)
        integer(kind=8) :: nnz(1:neq)
        integer(kind=8) :: qsize(neq)
        integer(kind=8) :: llist(neq)
        integer(kind=8) :: suiv(neq)
        integer(kind=8) :: p(neq)
        integer(kind=8) :: q(n2)
        integer(kind=8) :: invp(neq)
        integer(kind=8) :: perm(neq)
        integer(kind=8) :: lgind
        integer(kind=8) :: ddlmoy
        integer(kind=8) :: nbsn
        integer(kind=8) :: optnum
        integer(kind=8) :: nrl
        integer(kind=8) :: deb(*)
        integer(kind=8) :: vois(*)
        integer(kind=8) :: suit(*)
        integer(kind=8) :: ier
        integer(kind=8) :: nec
        integer(kind=8) :: prno(*)
        integer(kind=8) :: deeq(*)
        integer(kind=8) :: noeud(*)
        integer(kind=8) :: ddl(*)
        integer(kind=8) :: invpnd(*)
        integer(kind=8) :: permnd(*)
        integer(kind=8) :: spndnd(*)
        integer(kind=8) :: xadjd(*)
        aster_logical :: matgen
    end subroutine preml1
end interface
