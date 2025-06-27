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
                    ff_co, ff_3d, s, t, n, epai, crig, & 
                    nno_co, nno_3d, skip, mat )
    use raco3d_module
    !
    real(kind=8), intent(in) :: epai, crig
    integer(kind=8), intent(in) :: nb_gauss, nno_co, nno_3d
    real(kind=8), intent(in) :: jac_det(NB_GAUSS_MAX)
    real(kind=8), intent(in) :: gauss_weight(NB_GAUSS_MAX)
    real(kind=8), intent(in) :: gauss_coor(2, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: ff_co(NB_NO_CO_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: ff_3d(NB_NO_3D_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: t(3, NB_GAUSS_MAX), n(3, NB_GAUSS_MAX), s(3)
    aster_logical, intent(in) :: skip(NB_GAUSS_MAX)
    real(kind=8), intent(out) :: mat(:,:)
    end subroutine rco3d_calcmat
end interface
