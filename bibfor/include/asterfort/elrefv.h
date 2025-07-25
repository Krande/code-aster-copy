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
    subroutine elrefv(fami    , ndim    ,&
                      nnoL    , nnoQ    , nnos,&
                      npg     , jv_poids,&
                      jv_vfL  , jv_vfQ  ,&
                      jv_dfdeL, jv_dfdeQ,&
                      jv_ganoL, jv_ganoQ)
        character(len=4), intent(in) :: fami
        integer(kind=8), intent(out) :: ndim, nnos
        integer(kind=8), intent(out) :: npg, jv_poids
        integer(kind=8), intent(out) :: nnoL, jv_vfL, jv_dfdeL, jv_ganoL
        integer(kind=8), intent(out) :: nnoQ, jv_vfQ, jv_dfdeQ, jv_ganoQ
    end subroutine elrefv
end interface
