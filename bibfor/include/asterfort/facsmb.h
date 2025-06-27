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
    subroutine facsmb(nbnd, nbsn, supnd, invsup, parent,&
                      xadj, adjncy, anc, nouv, fils,&
                      frere, local, global, adress, lfront,&
                      nblign, lgsn, debfac, debfsn, chaine,&
                      place, nbass, delg, lgind, ier)
        integer(kind=8) :: lgind
        integer(kind=8) :: nbsn
        integer(kind=8) :: nbnd
        integer(kind=8) :: supnd(nbsn+1)
        integer(kind=8) :: invsup(nbnd)
        integer(kind=8) :: parent(nbsn)
        integer(kind=8) :: xadj(nbnd+1)
        integer(kind=8) :: adjncy(*)
        integer(kind=8) :: anc(nbnd)
        integer(kind=8) :: nouv(nbnd)
        integer(kind=8) :: fils(nbsn)
        integer(kind=8) :: frere(nbsn)
        integer(kind=4) :: local(lgind)
        integer(kind=4) :: global(lgind)
        integer(kind=8) :: adress(nbsn+1)
        integer(kind=8) :: lfront(nbsn)
        integer(kind=8) :: nblign(nbsn)
        integer(kind=8) :: lgsn(nbsn)
        integer(kind=8) :: debfac(nbnd+1)
        integer(kind=8) :: debfsn(nbsn+1)
        integer(kind=8) :: chaine(nbnd)
        integer(kind=8) :: place(nbnd)
        integer(kind=8) :: nbass(nbsn)
        integer(kind=8) :: delg(nbnd)
        integer(kind=8) :: ier
    end subroutine facsmb
end interface
