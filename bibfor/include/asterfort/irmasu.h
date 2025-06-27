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
    subroutine irmasu(ifc, ndim, nno, coordo, nbma,&
                      connex, point, typma, typel, codgra,&
                      codphy, codphd, permut, maxnod, lmod,&
                      noma, nbgrn, nogn, nbgrm, nogm,&
                      lmasu, nomai, nonoe, versio)
        integer(kind=8) :: maxnod
        integer(kind=8) :: ifc
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        real(kind=8) :: coordo(*)
        integer(kind=8) :: nbma
        integer(kind=8) :: connex(*)
        integer(kind=8) :: point(*)
        integer(kind=8) :: typma(*)
        integer(kind=8) :: typel(*)
        integer(kind=8) :: codgra(*)
        integer(kind=8) :: codphy(*)
        integer(kind=8) :: codphd(*)
        integer(kind=8) :: permut(maxnod, *)
        aster_logical :: lmod
        character(len=8) :: noma
        integer(kind=8) :: nbgrn
        character(len=24) :: nogn(*)
        integer(kind=8) :: nbgrm
        character(len=24) :: nogm(*)
        aster_logical :: lmasu
        character(len=8) :: nomai(*)
        character(len=8) :: nonoe(*)
        integer(kind=8) :: versio
    end subroutine irmasu
end interface
