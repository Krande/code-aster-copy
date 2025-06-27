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
    subroutine lckimp(ndim, typmod, option, mat, eps,&
                  phitot, vim, sig, forc_endo, vip, dsde_1, dsde_2, dsde_3)

    integer(kind=8)           :: ndim
    character(len=8)  :: typmod
    character(len=16) :: option
    integer(kind=8)           :: mat
    real(kind=8)      :: eps(:)
    real(kind=8)      :: phitot
    real(kind=8)      :: vim(:)
    real(kind=8)      :: vip(:)
    real(kind=8)      :: sig(:)
    real(kind=8)      :: forc_endo
    real(kind=8)      :: dsde_1(:,:)
    real(kind=8)      :: dsde_2(:)
    real(kind=8)      :: dsde_3

    end subroutine lckimp
end interface
