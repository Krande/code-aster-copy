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
#include "asterf_types.h"
interface 
      subroutine hydrxmat(xmat0, xmat1, hydra1, hydras, n,&
                          erreur)
        real(kind=8), intent(in) :: xmat0
        real(kind=8), intent(out) :: xmat1
        real(kind=8), intent(in) :: hydra1
        real(kind=8), intent(in) :: hydras
        real(kind=8), intent(in) :: n
        integer(kind=8), intent(out) :: erreur
    end subroutine hydrxmat
end interface
