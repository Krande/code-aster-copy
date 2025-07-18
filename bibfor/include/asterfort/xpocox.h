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
    subroutine xpocox(nbmac, ima, inmtot, nbcmpc, jresd1,&
                      jresv1, jresl1, jresd2, jresv2, jresl2)
        integer(kind=8) :: nbmac
        integer(kind=8) :: ima
        integer(kind=8) :: inmtot
        integer(kind=8) :: nbcmpc
        integer(kind=8) :: jresd1
        integer(kind=8) :: jresv1
        integer(kind=8) :: jresl1
        integer(kind=8) :: jresd2
        integer(kind=8) :: jresv2
        integer(kind=8) :: jresl2
    end subroutine xpocox
end interface
