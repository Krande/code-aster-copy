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
subroutine latrco(iTria, nbPoinInte, poinInte, triaCoorPara)
!
    implicit none
!
#include "asterfort/mesh_pairing_type.h"
!
    integer(kind=8), intent(in) :: iTria
    integer(kind=8), intent(in) :: nbPoinInte
    real(kind=8), intent(in) :: poinInte(2, MAX_NB_INTE)
    real(kind=8), intent(out) :: triaCoorPara(2, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (LAC) - Elementary computations
!
! Coordinates of current triangle
!
! --------------------------------------------------------------------------------------------------
!
! In  iTria            : index of current triangle
! In  nbPoinInte       : number of intersection points
! In  poinInte         : coordinates of intersection points (in cell parametric space)
! Out triaCoorPara     : coordinates of current triangle (in cell parametric space)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iPoinInte, iPoinInte1, iPoinInte2
    real(kind=8) ::  barycenter(2)
!
! --------------------------------------------------------------------------------------------------
!
    if (nbPoinInte == 3) then
        triaCoorPara(1:2, 1:3) = poinInte(1:2, 1:3)
    else
        barycenter = 0.d0
        do iPoinInte = 1, nbPoinInte
            barycenter(1:2) = barycenter(1:2)+poinInte(1:2, iPoinInte)
        end do
        barycenter = barycenter/real(nbPoinInte, kind=8)
        iPoinInte1 = iTria
        iPoinInte2 = iTria+1
        if (iTria == nbPoinInte) then
            iPoinInte2 = 1
        end if
        triaCoorPara(1:2, 1) = poinInte(1:2, iPoinInte1)
        triaCoorPara(1:2, 2) = poinInte(1:2, iPoinInte2)
        triaCoorPara(1:2, 3) = barycenter
    end if
!
end subroutine
