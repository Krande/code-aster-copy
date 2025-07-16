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
subroutine te0318(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_eval_module
!
    implicit none
#include "asterf_types.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! CALCUL DU FLUX AU CARRE AUX POINTS DE GAUSS
! ELEMENTS ISOPARAMETRIQUES OPTION : 'SOUR_ELGA '
!
!
! IN  OPTION : OPTION DE CALCUL
! IN  NOMTE  : NOM DU TYPE ELEMENT
!
! ----------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: nomres(1)
    character(len=32) :: phenom
    real(kind=8) :: lambda
    real(kind=8) :: valres(1), dtpg(3), dtpg_moy(3), dtpg_norm
    integer(kind=8) :: kp, itemps, iflux, imate
    real(kind=8), pointer :: tempe(:) => null()
! ------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPER', 'L', vr=tempe)
    call jevech('PSOUR_R', 'E', iflux)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
!
    if (phenom .eq. 'THER') then
        nomres(1) = 'LAMBDA'
    else
        call utmess('F', 'ELEMENTS2_67')
    end if
!
    dtpg_moy = 0.d0
!
    do kp = 1, FEQuadCell%nbQuadPoints
        call rcvalb('RIGI', kp, 1, '+', zi(imate), &
                    ' ', phenom, 1, 'INST', [zr(itemps)], &
                    1, nomres, valres, icodre, 1)
        lambda = valres(1)
        dtpg = FEEvalGradVec(FEBasis, tempe, FEQuadCell%points_param(1:3, kp))
        dtpg_moy = dtpg_moy-lambda*dtpg
    end do
!
    dtpg_moy = dtpg_moy/FEQuadCell%nbQuadPoints
    dtpg_norm = (dtpg_moy(1)**2+dtpg_moy(2)**2+dtpg_moy(3)**2)
    do kp = 1, FEQuadCell%nbQuadPoints
        call rcvalb('RIGI', kp, 1, '+', zi(imate), &
                    ' ', phenom, 1, 'INST', [zr(itemps)], &
                    1, nomres, valres, icodre, 1)
        zr(iflux+(kp-1)) = dtpg_norm/lambda
    end do
end subroutine
