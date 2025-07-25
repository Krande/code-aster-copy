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

    subroutine xsolveurtria(coor_nod, x, y, z, D, indmax, solution )
    
        integer(kind=8)                           ::  indmax    
        real(kind=8),dimension(3,3)       ::  coor_nod
        real(kind=8)                      ::  D(:)
        real(kind=8)                      ::  x(:)
        real(kind=8)                      ::  y(:)
        real(kind=8)                      ::  z(:)
        real(kind=8)                      ::  solution                                
   end subroutine
   
end interface
