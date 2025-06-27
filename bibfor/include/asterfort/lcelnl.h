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
    subroutine lcelnl(BEHinteg,&
                  fami, kpg, ksp, ndim, &
                  typmod, imate, compor, crit,&
                  option, eps, sig, vi, dsidep, codret)

    use Behaviour_type

    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*), option
    integer(kind=8) :: kpg, ksp, ndim, imate, codret
    real(kind=8) :: crit(*)
    real(kind=8) :: eps(:), sig(:), vi(1), dsidep(:,:)

    end subroutine lcelnl
end interface
