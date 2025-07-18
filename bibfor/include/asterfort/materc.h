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
    subroutine materc(matmas, matrig, matamo, numnu, amor, nommes,&
                      lfreqs, nbfreq,matobs, obsdim, gamma, alpha,eval)
#include "asterf_types.h"       
        character(len=8),intent(out) :: matmas
        character(len=8),intent(out) :: matrig
        character(len=8),intent(out) :: matamo
        character(len=8),intent(out) :: numnu
        aster_logical,intent(out) :: amor
        character(len=8),intent(out) :: nommes
        character(len=24),intent(out) :: lfreqs
        integer(kind=8),intent(out) :: nbfreq
        character(len=24),intent(out) :: matobs(3)
        integer(kind=8),intent(out) :: obsdim(3)
        real(kind=8),intent(out) :: gamma, alpha
        aster_logical,intent(out) :: eval
    end subroutine materc
end interface 
