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
subroutine xcalculfmm(nbno, jcalculs, jcopiels, jnodto, ndim, nodvois, &
                      jltno, jvcn, jgrlr, jbl, jbeta, jlistp , jvp, vale, &
                      deltat, levset, signls)
        integer(kind=8)           :: nbno                     
        integer(kind=8)           :: jcalculs
        integer(kind=8)           :: jcopiels
        integer(kind=8)           :: jnodto
        integer(kind=8)           :: ndim
        integer(kind=8)           :: nodvois 
        integer(kind=8)           :: jltno
        integer(kind=8)           :: jvcn
        integer(kind=8)           :: jgrlr
        integer(kind=8)           :: jbl
        integer(kind=8)           :: jbeta
        integer(kind=8)           :: jlistp
        integer(kind=8)           :: jvp
        real(kind=8)      :: vale(:)                                       
        real(kind=8)      :: deltat
        character(len=2)  :: levset
        character(len=3)  :: signls 
   end subroutine
end interface
