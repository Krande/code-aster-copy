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
    subroutine irmare(ifc, ndim, nno, coordo, nbma,&
                      connex, point, noma, typma, typel,&
                      lmod, titre, nbtitr, nbgrn, nbgrm,&
                      nomai, nonoe, formar)
        integer(kind=8) :: ifc
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        real(kind=8) :: coordo(*)
        integer(kind=8) :: nbma
        integer(kind=8) :: connex(*)
        integer(kind=8) :: point(*)
        character(len=8) :: noma
        integer(kind=8) :: typma(*)
        integer(kind=8) :: typel(*)
        aster_logical :: lmod
        character(len=80) :: titre(*)
        integer(kind=8) :: nbtitr
        integer(kind=8) :: nbgrn
        integer(kind=8) :: nbgrm
        character(len=8) :: nomai(*)
        character(len=8) :: nonoe(*)
        character(len=16) :: formar
    end subroutine irmare
end interface
