# Header only

In the linking phase of `bibfor` library on windows, I keep seeing these errors:

`error LNK2001: unresolved external symbol R8GAEM`

I've found this function declared in a fortran header file, but I can't seem to find the code itself. 
Is it possible that this is a header only fortran function?.
This compiles and links just fine using gfortran on linux.


```fortran
! --------------------------------------------------------------------
! Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org
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
function r8gaem()
real(kind=8) :: r8gaem
end function r8gaem
end interface
```