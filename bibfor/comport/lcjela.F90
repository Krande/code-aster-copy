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

subroutine lcjela(rela_comp, mod, nmat, mater, vin, &
                  dsde)
    implicit none
!       MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT ELASTIQUE A T+DT OU T
!       IN  LOI    :  MODELE DE COMPORTEMENT
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           VIN    :  VARIABLES INTERNES
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT
!       ----------------------------------------------------------------
#include "asterfort/lcopli.h"
#include "asterfort/rslpli.h"
    integer(kind=8) :: nmat
    real(kind=8) :: dsde(6, 6)
    real(kind=8) :: vin(*)
    real(kind=8) :: mater(nmat, 2)
    character(len=8) :: mod
    character(len=16), intent(in) :: rela_comp
!       ----------------------------------------------------------------
    if (rela_comp .eq. 'ROUSS_PR' .or. rela_comp .eq. 'ROUSS_VISC') then
        call rslpli('ISOTROPE', mod, mater, dsde, nmat, &
                    vin)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
!
        if (mater(nmat, 1) .eq. 0) then
            call lcopli('ISOTROPE', mod, mater(1, 1), dsde)
        else if (mater(nmat, 1) .eq. 1) then
            call lcopli('ORTHOTRO', mod, mater(1, 1), dsde)
        end if
!
!    CAS GENERAL : ELASTICITE LINEAIRE ISOTROPE OU ANISOTROPE
    else
!
        call lcopli('ISOTROPE', mod, mater(1, 1), dsde)
!
    end if
!
end subroutine
