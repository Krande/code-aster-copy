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

subroutine mm_cycl_d1_ss(pres_near, laug_cont_prev, laug_cont_curr, zone_cont_prev, &
                         zone_cont_curr, cycl_sub_type, alpha_cont_matr, alpha_cont_vect)
!
    implicit none
!
#include "asterfort/mm_cycl_zonc.h"
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
    real(kind=8), intent(in) :: pres_near
    real(kind=8), intent(in) :: laug_cont_prev
    real(kind=8), intent(in) :: laug_cont_curr
    real(kind=8), intent(out) :: alpha_cont_matr, alpha_cont_vect
    integer(kind=8), intent(out) :: zone_cont_prev
    integer(kind=8), intent(out) :: zone_cont_curr
    integer(kind=8), intent(out) :: cycl_sub_type
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve - Cycling
!
! Detection: contact/no-contact sub-cycling
!
! --------------------------------------------------------------------------------------------------
!
! In  pres_near        : tolerance for "near" contact - pressure
! In  laug_cont_prev   : previous augmented lagrangien for contact
! In  laug_cont_curr   : current augmented lagrangien for contact
! Out zone_cont_prev   : previous zone of contact
! Out zone_cont_curr   : current zone of contact
! Out cycl_sub_type    : sub-cycling type
!
! --------------------------------------------------------------------------------------------------
!
    cycl_sub_type = 0
    zone_cont_prev = 0
    zone_cont_curr = 0
!
! - Zoning detection - Previous
!
    call mm_cycl_zonc(pres_near, laug_cont_prev, zone_cont_prev)
!
! - Zoning detection - Current
!
    call mm_cycl_zonc(pres_near, laug_cont_curr, zone_cont_curr)
!
! - Sub-cycling 1 : grazing cycling
!
!   alpha_cont_matr = 0.3
!   alpha_cont_vect = 0.9
    if (((zone_cont_prev .eq. 3) .and. (zone_cont_curr .eq. 2)) .or. &
        ((zone_cont_prev .eq. 2) .and. (zone_cont_curr .eq. 3))) then
        cycl_sub_type = 1
!        alpha_cont_matr = 0.5
!        alpha_cont_vect = 1.0
        if (zone_cont_prev .eq. 3) then
            alpha_cont_matr = 0.95
            alpha_cont_vect = 0.95
        else
            alpha_cont_matr = 0.05
            alpha_cont_vect = 0.05
        end if

!
! - Sub-cycling 2
!
    elseif (((zone_cont_prev .eq. 2) .and. (zone_cont_curr .eq. 4)) .or. &
            ((zone_cont_prev .eq. 4) .and. (zone_cont_curr .eq. 2))) then
        cycl_sub_type = 2
        if (zone_cont_prev .eq. 4) then
            alpha_cont_matr = 1.0
            alpha_cont_vect = 1.0
        else
            alpha_cont_matr = 0.05
            alpha_cont_vect = 0.05
        end if

!
! - Sub-cycling 3
!
    elseif (((zone_cont_prev .eq. 1) .and. (zone_cont_curr .eq. 3)) .or. &
            ((zone_cont_prev .eq. 3) .and. (zone_cont_curr .eq. 1))) then
        cycl_sub_type = 3
        alpha_cont_matr = 1.0
        alpha_cont_vect = 1.0

!
! - Sub-cycling 4
!
    elseif (((zone_cont_prev .eq. 1) .and. (zone_cont_curr .eq. 4)) .or. &
            ((zone_cont_prev .eq. 4) .and. (zone_cont_curr .eq. 1))) then
        cycl_sub_type = 4
        alpha_cont_matr = 0.95
        alpha_cont_vect = 1.0

    else
        cycl_sub_type = 5
        alpha_cont_matr = 1.0
        alpha_cont_vect = 1.0
    end if

end subroutine
