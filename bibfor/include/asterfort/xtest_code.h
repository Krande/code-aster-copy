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
interface
    function xtest_code(id1, id2, lfno, nfh, nfissmax, fno1, fno2)
       integer(kind=8) :: id1
       integer(kind=8) :: id2
       integer(kind=8) :: nfh
       integer(kind=8) :: nfissmax
       aster_logical :: xtest_code
       aster_logical :: lfno
       integer(kind=8) :: fno1(nfissmax)
       integer(kind=8) :: fno2(nfissmax)
    end function xtest_code
end interface
