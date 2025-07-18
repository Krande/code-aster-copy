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
    subroutine rcfonc(quest, ktrac , jprol , jvale, nbvale,&
                      sigy , e     , nu    , p    , rp    ,&
                      rprim, airerp, sieleq, dp)
        character(len=1), intent(in) :: quest
        integer(kind=8), intent(in) :: ktrac
        integer(kind=8), intent(in) :: jprol
        integer(kind=8), intent(in) :: jvale
        integer(kind=8), intent(in) :: nbvale
        real(kind=8), optional, intent(in) :: e
        real(kind=8), optional, intent(in) :: nu
        real(kind=8), optional, intent(in) :: sieleq
        real(kind=8), optional, intent(in) :: p
        real(kind=8), optional, intent(out) :: sigy
        real(kind=8), optional, intent(out) :: rp
        real(kind=8), optional, intent(out) :: rprim
        real(kind=8), optional, intent(out) :: airerp
        real(kind=8), optional, intent(out) :: dp
    end subroutine rcfonc
end interface
