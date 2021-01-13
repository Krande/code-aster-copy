! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
    subroutine xtform(typmae, typmam, typmac,&
                      nnm, coore, coorm, coorc,&
                      ffe, ffm, dffc)
        character(len=8), intent(in) :: typmae, typmam, typmac
        real(kind=8), intent(in) :: coorc(2), coore(3), coorm(3)
        integer, intent(in) :: nnm
        real(kind=8), intent(out) :: ffe(20), ffm(20), dffc(3, 9)
    end subroutine xtform
end interface
