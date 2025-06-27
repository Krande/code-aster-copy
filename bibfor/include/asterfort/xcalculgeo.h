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
   subroutine xcalculgeo(ndim, vale, jvp, jbl, deltat, jbeta, &
                         jlistp, node, newlst, newlsn)
       integer(kind=8)          :: ndim                  
       real(kind=8)     :: vale(:)
       integer(kind=8)          :: jvp
       integer(kind=8)          :: jbl
       real(kind=8)     :: deltat
       integer(kind=8)          :: jbeta
       integer(kind=8)          :: jlistp
       integer(kind=8)          :: node                                
       real(kind=8)     :: newlst 
       real(kind=8)     :: newlsn
   end subroutine
end interface
