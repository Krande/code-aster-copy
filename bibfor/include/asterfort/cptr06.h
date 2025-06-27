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
    subroutine cptr06(maout , inc   , jcnnpa, conloc,&
                      limane, jmacou, jmacsu, macou ,&
                      macsu , ind   , ind1  )
        character(len=8), intent(in) :: maout
        integer(kind=8), intent(in) :: inc 
        integer(kind=8), intent(in) :: jcnnpa
        character(len=24), intent(in) :: conloc
        character(len=24), intent(in) :: limane
        integer(kind=8), intent(in) :: jmacou
        integer(kind=8), intent(in) :: jmacsu
        integer(kind=8), intent(in) :: macou
        integer(kind=8), intent(in) :: macsu
        integer(kind=8), intent(out) :: ind 
        integer(kind=8), intent(out) :: ind1 
    end subroutine cptr06
end interface
