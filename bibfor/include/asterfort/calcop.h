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
    subroutine calcop(option, listOptJvZ, resultIn, resultOut, listStoreJv, &
                      nbStore, resultType, codret, jvBase_, tldist_)
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: listOptJvZ
        character(len=8), intent(in) :: resultIn, resultOut
        character(len=19), intent(in) :: listStoreJv
        integer(kind=8), intent(in) :: nbStore
        character(len=16), intent(in) :: resultType
        integer(kind=8), intent(out) ::  codret
        character(len=1), optional, intent(in) :: jvBase_
        aster_logical, optional, intent(in)  :: tldist_
    end subroutine calcop
end interface
