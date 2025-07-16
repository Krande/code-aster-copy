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
#include "MeshTypes_type.h"
!
interface
    subroutine irmmma(fid, nomamd, nbCell, connex, point,&
                      typma, nommai, prefix, nbtyp, typgeo,&
                      nomtyp, nnotyp, renumd, nbCellType, infmed,&
                      modnum, nuanom, nosdfu)
        med_idt :: fid
        integer(kind=8) :: nbCell, nbtyp
        integer(kind=8) :: connex(*), typma(*), point(*)
        integer(kind=8) :: typgeo(*), nnotyp(*), nbCellType(MT_NTYMAX)
        integer(kind=8) :: renumd(*), modnum(MT_NTYMAX), nuanom(MT_NTYMAX, *)
        integer(kind=8) :: infmed
        character(len=6) :: prefix
        character(len=8) :: nommai(*)
        character(len=8) :: nomtyp(*)
        character(len=*) :: nomamd
        character(len=8) :: nosdfu
    end subroutine irmmma
end interface
