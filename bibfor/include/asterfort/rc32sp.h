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
    subroutine rc32sp(ze200, lieu, iocc1, iocc2, ns,&
                      sp, spmeca, instsp, nbsscyc, spss)
        aster_logical :: ze200
        character(len=4) :: lieu
        integer(kind=8) :: iocc1
        integer(kind=8) :: iocc2
        integer(kind=8) :: ns
        real(kind=8) :: sp(2)
        real(kind=8) :: spmeca(2)
        real(kind=8) :: instsp(4)
        integer(kind=8) :: nbsscyc
        real(kind=8) :: spss(100)
    end subroutine rc32sp
end interface
