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
#include "asterf_types.h"
!
interface
    subroutine mmmvcc(phase , l_pena_cont, &
                      nnl   , wpg        , ffl   , jacobi,&
                      jeu   , coefac     , dlagrc,&
                      vectcc)
        character(len=4), intent(in) :: phase
        aster_logical, intent(in) :: l_pena_cont
        integer(kind=8), intent(in) :: nnl
        real(kind=8), intent(in) :: wpg, ffl(9), jacobi
        real(kind=8), intent(in) :: jeu, dlagrc, coefac
        real(kind=8), intent(out) :: vectcc(9)
    end subroutine mmmvcc
end interface
