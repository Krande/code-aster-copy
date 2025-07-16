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

subroutine te0078(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_rhs_module
    use FE_eval_module
!
    implicit none
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/nlcomp.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
!
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_EVOL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
!
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadRigi, FEQuadMass
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: nbres
    parameter(nbres=1)
    integer(kind=8) :: icodre(nbres)
    character(len=16) :: phenom
    real(kind=8) :: valQPM(MAX_QP), tpg, dtpg(3), flux(3), BGSEval(3, MAX_BS)
    real(kind=8) :: resi_f(MAX_BS), resi_m(MAX_BS), resi(MAX_BS)
    real(kind=8) :: cp, valres(1), Kglo(3, 3), time, deltat, theta
    integer(kind=8) :: kp, imate, icamas, itemps
    real(kind=8), pointer :: temp(:) => null()
    character(len=8), parameter :: famiR = "RIGI"
    character(len=8), parameter :: famiM = "MASS"
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    call FEQuadMass%initCell(FECell, famiM)
    call FEQuadRigi%initCell(FECell, famiR)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPER', 'L', vr=temp)
!
    time = zr(itemps)
    deltat = zr(itemps+1)
    theta = zr(itemps+2)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
!
!   pour stopper le calcul si PCAMASS n'est pas disponible
    if (phenom == "THER_ORTH") then
        call jevech('PCAMASS', 'L', icamas)
    end if
!
    resi_f = 0.d0
    do kp = 1, FEQuadRigi%nbQuadPoints
        BGSEval = FEBasis%grad(FEQuadRigi%points_param(1:3, kp), FEQuadRigi%jacob(1:3, 1:3, kp))
!
        dtpg = FEEvalGradVec(FEBasis, temp, FEQuadRigi%points_param(1:3, kp), BGSEval)
!
        call nlcomp(phenom, famiR, kp, imate, FECell%ndim, FEQuadRigi%points(1:3, kp), &
                    time, 0.d0, Kglo, dtp_=dtpg, fluglo_=flux)
!
        call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadRigi%weights(kp), flux, resi_f)
    end do
!
    do kp = 1, FEQuadMass%nbQuadPoints
!
        call rcvalb(famiM, kp, 1, '+', zi(imate), ' ', phenom, 1, 'INST', [time], &
                    1, 'RHO_CP', valres, icodre, 1)
        cp = valres(1)
!
        tpg = FEEvalFuncRScal(FEBasis, temp, FEQuadMass%points_param(1:3, kp))
        ValQPM(kp) = cp*tpg
    end do
!
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPM, resi_m)
!
    resi = (theta-1.0d0)*resi_f+resi_m/deltat
!
    call writeVector('PVECTTR', FEBasis%size, resi)
!
end subroutine
