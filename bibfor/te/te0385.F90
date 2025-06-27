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

subroutine te0385(nomopt, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_eval_module
!
    implicit none
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/runge6.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: nomte, nomopt
! ----------------------------------------------------------------------
!
!    - FONCTION REALISEE:  OPTION : 'HYDR_ELGA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! THERMIQUE NON LINEAIRE
!
! ----------------------------------------------------------------------
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: itemps, imate
    integer(kind=8) ::  ifon(6), kp
    real(kind=8) :: deltat, err, tpgm, tpgp
    real(kind=8), pointer :: hydrgm(:) => null(), hydrgp(:) => null()
    real(kind=8), pointer :: tempm(:) => null()
    real(kind=8), pointer :: tempp(:) => null()
    character(len=16), pointer :: compor(:) => null()
!
    if (nomopt .ne. "HYDR_ELGA") then
        ASSERT(ASTER_FALSE)
    end if
!
    call jevech('PCOMPOR', 'L', vk16=compor)

    if (compor(RELA_NAME) (1:9) .eq. 'THER_HYDR') then
        call FECell%init()
        call FEBasis%initCell(FECell)
        call FEQuadCell%initCell(FECell, "MASS")
!
        call jevech('PMATERC', 'L', imate)
        call jevech('PTEMPMR', 'L', vr=tempm)
        call jevech('PTEMPPR', 'L', vr=tempp)
        call jevech('PHYDRMR', 'L', vr=hydrgm)
        call jevech('PHYDRPR', 'E', vr=hydrgp)
        call jevech('PINSTR', 'L', itemps)
!
        deltat = zr(itemps+1)
!
        call ntfcma(compor(RELA_NAME), zi(imate), ASTER_FALSE, ifon)
!
        do kp = 1, FEQuadCell%nbQuadPoints
            hydrgp(kp) = 0.d0
            tpgm = FEEvalFuncRScal(FEBasis, tempm, FEQuadCell%points_param(1:3, kp))
            tpgp = FEEvalFuncRScal(FEBasis, tempp, FEQuadCell%points_param(1:3, kp))
!
            call runge6(ifon(3), deltat, tpgp, tpgm, hydrgm(kp), &
                        hydrgp(kp), err)
        end do
    end if
end subroutine
