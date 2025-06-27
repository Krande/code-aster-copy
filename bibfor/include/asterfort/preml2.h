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
    subroutine preml2(n1, diag, col, delg, xadj1,&
                      adjnc1, estim, adress, parend, fils,&
                      frere, anc, nouv, supnd, dhead,&
                      qsize, llist, marker, invsup, local,&
                      global, lfront, nblign, decal, lgsn,&
                      debfac, debfsn, seq, lmat, adpile,&
                      chaine, suiv, place, nbass, ncbloc,&
                      lgbloc, nbloc, lgind, nbsnd, ier)
        integer(kind=8) :: n1
        integer(kind=8) :: diag(0:n1)
        integer(kind=8) :: col(*)
        integer(kind=8) :: delg(*)
        integer(kind=8) :: xadj1(n1+1)
        integer(kind=8) :: adjnc1(*)
        integer(kind=8) :: estim
        integer(kind=8) :: adress(*)
        integer(kind=8) :: parend(*)
        integer(kind=8) :: fils(n1)
        integer(kind=8) :: frere(n1)
        integer(kind=8) :: anc(n1)
        integer(kind=8) :: nouv(n1)
        integer(kind=8) :: supnd(n1)
        integer(kind=8) :: dhead(*)
        integer(kind=8) :: qsize(*)
        integer(kind=8) :: llist(*)
        integer(kind=8) :: marker(*)
        integer(kind=8) :: invsup(n1)
        integer(kind=4) :: local(*)
        integer(kind=4) :: global(*)
        integer(kind=8) :: lfront(n1)
        integer(kind=8) :: nblign(n1)
        integer(kind=8) :: decal(*)
        integer(kind=8) :: lgsn(n1)
        integer(kind=8) :: debfac(n1+1)
        integer(kind=8) :: debfsn(n1)
        integer(kind=8) :: seq(n1)
        integer(kind=8) :: lmat
        integer(kind=8) :: adpile(n1)
        integer(kind=8) :: chaine(n1)
        integer(kind=8) :: suiv(n1)
        integer(kind=8) :: place(n1)
        integer(kind=8) :: nbass(n1)
        integer(kind=8) :: ncbloc(*)
        integer(kind=8) :: lgbloc(*)
        integer(kind=8) :: nbloc
        integer(kind=8) :: lgind
        integer(kind=8) :: nbsnd
        integer(kind=8) :: ier
    end subroutine preml2
end interface
