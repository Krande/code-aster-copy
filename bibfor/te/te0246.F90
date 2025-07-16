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
subroutine te0246(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_mass_module
    use FE_eval_module
!
    implicit none
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcfode.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_THER' ET 'MASS_THER_TANG'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: phenom, rela_name
    real(kind=8) :: valQP(MAX_QP), tpgi, r8bid
    real(kind=8) :: mass(MAX_BS, MAX_BS)
    integer(kind=8) :: kp, imate
    integer(kind=8) :: ifon(6)
    aster_logical :: aniso
    character(len=16), pointer :: compor(:) => null()
    real(kind=8), pointer :: tempi(:) => null()
!
!-----------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "MASS")
    call FEBasis%initCell(FECell)
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_name = compor(RELA_NAME)
    if (rela_name(1:5) .eq. 'THER_') then
        call jevech('PTEMPEI', 'L', vr=tempi)
        call jevech('PMATERC', 'L', imate)
!
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = ASTER_FALSE
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = ASTER_TRUE
        end if
        call ntfcma(rela_name, zi(imate), aniso, ifon)
    end if
!
    valQP = 0.0
    do kp = 1, FEQuadCell%nbQuadPoints
        if (rela_name(1:5) .eq. 'THER_') then
            tpgi = FEEvalFuncRScal(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
            call rcfode(ifon(1), tpgi, r8bid, valQP(kp))
        else if (rela_name(1:5) .eq. 'SECH_') then
            valQP(kp) = 1.d0
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call FEMassMatScal(FEQuadCell, FEBasis, mass, valQP)
!
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, mass)
!
end subroutine
