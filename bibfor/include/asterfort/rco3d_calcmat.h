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
subroutine rco3d_calcmat(nb_gauss, gauss_weight, gauss_coor, jac_det, &
                    ff_co, ff_3d, s, t, n, epai, & 
                    nno_co, nno_3d, mat )
    real(kind=8), intent(in) :: epai
    integer, intent(in) :: nb_gauss, nno_co, nno_3d
    real(kind=8), intent(in) :: jac_det(10)
    real(kind=8), intent(in) :: gauss_weight(10)
    real(kind=8), intent(in) :: gauss_coor(2, 10)
    real(kind=8), intent(in) :: ff_co(3, 10)
    real(kind=8), intent(in) :: ff_3d(8, 10)
    real(kind=8), intent(in) :: t(3, 10), n(3, 10), s(3)
    real(kind=8), intent(out) :: mat(:,:)
    end subroutine rco3d_calcmat
end interface

