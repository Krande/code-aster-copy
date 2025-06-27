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
    subroutine nmelru(fami, kpg, ksp, imate, compor, &
                      epseq, p_arg, divu, gonf, inco, nonlin, ener)
    
    character(len=*),intent(in)  :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    integer(kind=8), intent(in)          :: imate
    character(len=16),intent(in) :: compor(*)
    real(kind=8), intent(in)     :: epseq
    real(kind=8), intent(in)     :: p_arg
    real(kind=8), intent(in)     :: divu
    real(kind=8), intent(in)     :: gonf
    aster_logical, intent(in)    :: inco
    aster_logical, intent(in)    :: nonlin
    real(kind=8), intent(out)    :: ener(2)
    
    end subroutine 
end interface
