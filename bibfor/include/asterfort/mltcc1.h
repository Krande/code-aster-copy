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
    subroutine mltcc1(nbloc, ncbloc, decal, supnd, fils,&
                      frere, seq, lgsn, lfront, adress,&
                      local, adpile, nbass, pile, lgpile,&
                      adper, t1, t2, factol, factou,&
                      typsym, ad, eps, ier, nbb,&
                      cl, cu)
        integer(kind=8) :: nbb
        integer(kind=8) :: nbloc
        integer(kind=8) :: ncbloc(*)
        integer(kind=8) :: decal(*)
        integer(kind=8) :: supnd(*)
        integer(kind=8) :: fils(*)
        integer(kind=8) :: frere(*)
        integer(kind=8) :: seq(*)
        integer(kind=8) :: lgsn(*)
        integer(kind=8) :: lfront(*)
        integer(kind=8) :: adress(*)
        integer(kind=4) :: local(*)
        integer(kind=8) :: adpile(*)
        integer(kind=8) :: nbass(*)
        complex(kind=8) :: pile(*)
        integer(kind=8) :: lgpile
        integer(kind=8) :: adper(*)
        complex(kind=8) :: t1(*)
        complex(kind=8) :: t2(*)
        character(len=24) :: factol
        character(len=24) :: factou
        integer(kind=8) :: typsym
        integer(kind=8) :: ad(*)
        real(kind=8) :: eps
        integer(kind=8) :: ier
        complex(kind=8) :: cl(nbb, nbb, *)
        complex(kind=8) :: cu(nbb, nbb, *)
    end subroutine mltcc1
end interface
