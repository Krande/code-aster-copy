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

module searchlist_module

    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"

! -------------------------------------------------------
!
! Searchlist utilities
!
! -------------------------------------------------------

contains

    function getBounds(v, prec, crit) result(bounds)
        implicit none

        real(kind=8), intent(in) :: v, prec
        character(len=8), intent(in) :: crit
        real(kind=8) :: b1, b2, bounds(2)

        if (crit(1:7) .eq. 'RELATIF') then
            b1 = v*(1.d0-prec)
            b2 = v*(1.d0+prec)
        else if (crit(1:6) .eq. 'ABSOLU') then
            b1 = v-prec
            b2 = v+prec
        else
            ASSERT(ASTER_FALSE)
        end if

        bounds(1) = min(b1, b2)
        bounds(2) = max(b1, b2)

    end function getBounds

    function almostEqual(a, b, prec, crit) result(out)
        implicit none

        real(kind=8), intent(in) :: a, b, prec
        character(len=8), intent(in) :: crit
        aster_logical :: out
        !
        real(kind=8) :: bounds(2)
        !
        bounds = getBounds(b, prec, crit)
        out = (a >= bounds(1) .and. a <= bounds(2))

    end function almostEqual

    function searchCandidates(v, array, prec, crit, mask) result(out)
        implicit none

        real(kind=8), intent(in) :: v, prec, array(:)
        character(len=8), intent(in) :: crit
        aster_logical, intent(in) :: mask(size(array))
        integer(kind=8), allocatable :: out(:)
        !
        real(kind=8) :: bounds(2)
        integer(kind=8) :: i
        !
        bounds = getBounds(v, prec, crit)
        out = pack([(i, i=1, size(array))], &
                   (mask .and. (array >= bounds(1) .and. array <= bounds(2))))

    end function searchCandidates

    function isUnique(v, array, prec, crit, umask) result(out)
        implicit none

        real(kind=8), intent(in) :: v, prec, array(:)
        character(len=8), intent(in) :: crit
        aster_logical, optional, intent(in) :: umask(size(array))
        aster_logical :: out, mask(size(array))
        !
        if (.not. present(umask)) then
            mask = ASTER_TRUE
        else
            mask = umask
        end if
        !
        out = 1 .eq. size(searchCandidates(v, array, prec, crit, mask))

    end function isUnique

    function getUnique(v, array, prec, crit, umask) result(out)
        implicit none

        real(kind=8), intent(in) :: v, prec, array(:)
        character(len=8), intent(in) :: crit
        aster_logical, optional, intent(in) :: umask(size(array))
        aster_logical :: mask(size(array))
        integer(kind=8) :: out
        !
        integer(kind=8), allocatable :: candidates(:)
        !
        if (.not. present(umask)) then
            mask = ASTER_TRUE
        else
            mask = umask
        end if
        !
        out = 0
        candidates = searchCandidates(v, array, prec, crit, mask)
        if (size(candidates) .eq. 1) then
            out = candidates(1)
        end if

    end function getUnique

    function allUnique(array, prec, crit, umask) result(out)
        implicit none
        real(kind=8), intent(in) :: prec, array(:)
        character(len=8), intent(in) :: crit
        aster_logical, optional, intent(in) :: umask(size(array))
        aster_logical :: out, mask(size(array))
        !
        integer(kind=8) :: i
        !
        if (.not. present(umask)) then
            mask = ASTER_TRUE
        else
            mask = umask
        end if

        !
        out = ASTER_TRUE
        do i = 1, size(array)
            out = out .and. isUnique(array(i), array, prec, crit, mask)
        end do

    end function allUnique

end module searchlist_module
