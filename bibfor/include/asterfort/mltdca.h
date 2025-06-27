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
    subroutine mltdca(nbloc, lgbloc, ncbloc, decal, seq,&
                      nbsn, nbnd, supnd, adress, global,&
                      lgsn, factol, factou, sm, x,&
                      invp, perm, ad, trav, typsym)
        integer(kind=8) :: nbnd
        integer(kind=8) :: nbsn
        integer(kind=8) :: nbloc
        integer(kind=8) :: lgbloc(nbsn)
        integer(kind=8) :: ncbloc(nbnd)
        integer(kind=8) :: decal(nbsn)
        integer(kind=8) :: seq(nbsn)
        integer(kind=8) :: supnd(nbsn+1)
        integer(kind=8) :: adress(nbsn+1)
        integer(kind=4) :: global(*)
        integer(kind=8) :: lgsn(nbsn)
        character(len=24) :: factol
        character(len=24) :: factou
        complex(kind=8) :: sm(nbnd)
        complex(kind=8) :: x(nbnd)
        integer(kind=8) :: invp(nbnd)
        integer(kind=8) :: perm(nbnd)
        integer(kind=8) :: ad(nbnd)
        complex(kind=8) :: trav(nbnd)
        integer(kind=8) :: typsym
    end subroutine mltdca
end interface
