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
#include "asterf_types.h"
!
interface
    subroutine rco3d_apco3d(noma, lismavo, lismaco, nbmavo, nbmaco, epai, &
                        list_pairs, nb_pairs, nt_nodes)
        character(len=8), intent(in) :: noma
        character(len=24), intent(in) :: lismaco, lismavo
        integer(kind=8), intent(in) :: nbmavo, nbmaco
        real(kind=8), intent(in) :: epai
        integer(kind=8), intent(out) :: nb_pairs, nt_nodes
        integer(kind=8), pointer :: list_pairs(:)

    end subroutine rco3d_apco3d
end interface
