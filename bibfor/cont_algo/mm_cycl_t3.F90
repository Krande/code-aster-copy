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

subroutine mm_cycl_t3(pres_frot_prev, dist_frot_prev, coef_frot_prev, &
                      cycl_stat_curr)
!
    implicit none
!
!
! person_in_charge: mickael.abbas at edf.fr
!
    real(kind=8), intent(in) :: pres_frot_prev(3)
    real(kind=8), intent(in) :: dist_frot_prev(3)
    real(kind=8), intent(in) :: coef_frot_prev
    integer(kind=8), intent(out) :: cycl_stat_curr
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve - Cycling
!
! Treatment: sliding forward/backward
!
! --------------------------------------------------------------------------------------------------
!
! In  pres_frot_prev   : previous friction pressure in cycle
! In  dist_frot_prev   : previous friction distance in cycle
! In  coef_frot_prev   : previous augmented ratio for friction
! In  pres_frot_curr   : current friction pressure
! In  dist_frot_curr   : current friction distance
! Out cycl_stat_curr   : state of treatment
!                      -10 : Failure during adaptation
!                      -02 : Cannot adapt
!                      -01 : has been adapted
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: cycl_type, idim
    real(kind=8) :: laug_frot_curr(3), laug_frot_prev(3)
    real(kind=8) :: nrese_curr, nrese_prev
!
! --------------------------------------------------------------------------------------------------
!

!
! - Initialisations
!
    cycl_type = 3
    nrese_curr = 0.d0
    nrese_prev = 0.d0
!
! - Augmented ratios
!
    laug_frot_prev(1) = pres_frot_prev(1)+coef_frot_prev*dist_frot_prev(1)
    laug_frot_prev(2) = pres_frot_prev(2)+coef_frot_prev*dist_frot_prev(2)
    laug_frot_prev(3) = pres_frot_prev(3)+coef_frot_prev*dist_frot_prev(3)
    do idim = 1, 3
        nrese_prev = laug_frot_prev(idim)*laug_frot_prev(idim)+nrese_prev
    end do
    nrese_prev = sqrt(nrese_prev)
!
    laug_frot_curr(1) = pres_frot_prev(1)+coef_frot_prev*dist_frot_prev(1)
    laug_frot_curr(2) = pres_frot_prev(2)+coef_frot_prev*dist_frot_prev(2)
    laug_frot_curr(3) = pres_frot_prev(3)+coef_frot_prev*dist_frot_prev(3)
    do idim = 1, 3
        nrese_curr = laug_frot_curr(idim)*laug_frot_curr(idim)+nrese_curr
    end do
    nrese_curr = sqrt(nrese_curr)
    cycl_stat_curr = -2

end subroutine
