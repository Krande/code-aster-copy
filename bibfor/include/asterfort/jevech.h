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
    subroutine jevech(nmparz, louez, itab, vl, vi,&
                  vi4, vr, vc, vk8, vk16,&
                  vk24, vk32, vk80)
#include "asterf_types.h"
        character(len=*), intent(in) :: nmparz
        character(len=*), intent(in) :: louez
        integer(kind=8), optional :: itab

        aster_logical,     pointer, optional :: vl(:)
        integer(kind=8),           pointer, optional :: vi(:)
        integer(kind=4),   pointer, optional :: vi4(:)
        real(kind=8),      pointer, optional :: vr(:)
        complex(kind=8),   pointer, optional :: vc(:)
        character(len=8),  pointer, optional :: vk8(:)
        character(len=16), pointer, optional :: vk16(:)
        character(len=24), pointer, optional :: vk24(:)
        character(len=32), pointer, optional :: vk32(:)
        character(len=80), pointer, optional :: vk80(:)
    end subroutine jevech
end interface
