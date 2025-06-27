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
    subroutine mltpas(nbnd, nbsn, supnd, xadj, adjncy,&
                      anc, nouv, seq, global, adress,&
                      nblign, lgsn, nbloc, ncbloc, lgbloc,&
                      diag, col, lmat, place)
        integer(kind=8) :: nbsn
        integer(kind=8) :: nbnd
        integer(kind=8) :: supnd(nbsn+1)
        integer(kind=8) :: xadj(nbnd+1)
        integer(kind=8) :: adjncy(*)
        integer(kind=8) :: anc(nbnd)
        integer(kind=8) :: nouv(nbnd)
        integer(kind=8) :: seq(nbsn)
        integer(kind=4) :: global(*)
        integer(kind=8) :: adress(nbsn+1)
        integer(kind=8) :: nblign(nbsn)
        integer(kind=8) :: lgsn(nbsn)
        integer(kind=8) :: nbloc
        integer(kind=8) :: ncbloc(*)
        integer(kind=8) :: lgbloc(*)
        integer(kind=8) :: diag(0:nbnd)
        integer(kind=8) :: col(*)
        integer(kind=8) :: lmat
        integer(kind=8) :: place(nbnd)
    end subroutine mltpas
end interface
