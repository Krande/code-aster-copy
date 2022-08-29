! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine nmorth(fami, kpg, ksp, ndim, phenom,&
                      imate, poum, deps, sigm, option,&
                      angmas, sigp, dsidep)
                      
    character(len=*), intent(in) :: fami
    integer, intent(in)          :: kpg
    integer, intent(in)          :: ksp
    integer, intent(in)          :: ndim
    character(len=16), intent(in):: phenom
    integer, intent(in)          :: imate
    character(len=*), intent(in) :: poum
    real(kind=8),intent(in)      :: deps(2*ndim)
    real(kind=8),intent(in)      :: sigm(2*ndim)
    character(len=16), intent(in):: option
    real(kind=8),intent(in)      :: angmas(3)
    real(kind=8),intent(out)     :: sigp(2*ndim)
    real(kind=8),intent(out)     :: dsidep(2*ndim,2*ndim)
    
    end subroutine nmorth
end interface
