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
    subroutine rco3d_infos(typmaco, typma3d, epai, j_geom, nb_gauss, gauss_coor, &
                        gauss_weight, jac_det, ff_co, ff_3d, s, t, n, skip)
        use raco3d_module
        character(len=8), intent(in) :: typmaco, typma3d
        real(kind=8), intent(in) :: epai
        integer(kind=8), intent(in) :: j_geom
        integer(kind=8), intent(out) :: nb_gauss
        real(kind=8), intent(out) :: gauss_coor(2, NB_GAUSS_MAX)
        real(kind=8), intent(out) :: jac_det(NB_GAUSS_MAX)
        real(kind=8), intent(out) :: gauss_weight(NB_GAUSS_MAX)
        real(kind=8), intent(out) :: ff_co(NB_NO_CO_MAX,NB_GAUSS_MAX)
        real(kind=8), intent(out) :: ff_3d(NB_NO_3D_MAX, NB_GAUSS_MAX)
        real(kind=8), intent(out) :: t(3, NB_GAUSS_MAX), n(3, NB_GAUSS_MAX), s(3)
        aster_logical, intent(out) :: skip(NB_GAUSS_MAX)
    end subroutine rco3d_infos
end interface
