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
    subroutine lcvisc(fami, kpg, ksp, ndim, imate,&
                      lSigm, lMatr, lVari, &
                      instam, instap, deps, vim, &
                      sigp, vip, dsidep)
        character(len=*),intent(in) :: fami
        integer(kind=8),intent(in)          :: kpg
        integer(kind=8),intent(in)          :: ksp
        integer(kind=8),intent(in)          :: ndim
        integer(kind=8),intent(in)          :: imate
        aster_logical, intent(in)   :: lSigm
        aster_logical, intent(in)   :: lMatr
        aster_logical, intent(in)   :: lVari
        real(kind=8),intent(in)     :: instam
        real(kind=8),intent(in)     :: instap
        real(kind=8),intent(in)     :: deps(:)
        real(kind=8),intent(in)     :: vim(:)
        real(kind=8),intent(inout)  :: sigp(:)
        real(kind=8),intent(out)    :: vip(:)
        real(kind=8),intent(inout)  :: dsidep(:,:)
    end subroutine lcvisc
end interface
