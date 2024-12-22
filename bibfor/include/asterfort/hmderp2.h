! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
#include "asterf_types.h"
!
interface
    subroutine hmderp2(ds_thm, t, pvp, &
                      rho11 , rho12 , h11  , h12,&
                      dp11p1, dp11p2, dp11t,&
                      dp12p1, dp12p2, dp12t,&
                      dp21p1, dp21p2, dp21t,&
                      dp1pp1, dp2pp1, dtpp1,&
                      dp1pp2, dp2pp2, dtpp2,&
                      dp1pt , dp2pt , dtpt)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        real(kind=8), intent(in) :: t, pvp
        real(kind=8), intent(in) :: rho11, rho12, h11, h12
        real(kind=8), intent(out) :: dp11p1, dp11p2, dp11t
        real(kind=8), intent(out) :: dp12p1, dp12p2, dp12t
        real(kind=8), intent(out) :: dp21p1, dp21p2, dp21t
        real(kind=8), intent(out) :: dp1pp1, dp2pp1, dtpp1
        real(kind=8), intent(out) :: dp1pp2, dp2pp2, dtpp2
        real(kind=8), intent(out) :: dp1pt, dp2pt, dtpt
    end subroutine hmderp2
end interface
