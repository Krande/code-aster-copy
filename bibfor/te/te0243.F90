! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0243(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ntcomp.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'RAPH_THER'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer :: icamas, nbres
    parameter(nbres=3)
    integer :: icodre(nbres)
    character(len=32) :: phenom
    real(kind=8) ::  tpg, dtpg(3), tpsec, diff, fluglo(3)
    real(kind=8) :: resi(MAX_BS), valQP(3, MAX_QP)
    real(kind=8), pointer :: flux(:) => null()
    real(kind=8), pointer :: tempi(:) => null()
    real(kind=8), pointer :: sechf(:) => null()
    integer ::  kp, ifon(6)
    integer ::  imate, icomp, nbDof
    aster_logical :: aniso
! ----------------------------------------------------------------------
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
    nbDof = FEBasis%size
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPEI', 'L', vr=tempi)
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PFLUXPR', 'E', vr=flux)
!
    if ((zk16(icomp) (1:5) .eq. 'SECH_')) then
        if (zk16(icomp) (1:12) .eq. 'SECH_GRANGER' .or. zk16(icomp) (1:10) .eq. 'SECH_NAPPE') then
            call jevech('PTMPCHF', 'L', vr=sechf)
        else
!          POUR LES AUTRES LOIS, PAS DE CHAMP DE TEMPERATURE
!          ISECHF EST FICTIF
            call jevech('PTEMPEI', 'L', vr=sechf)
        end if
!
    else if (zk16(icomp) (1:5) .eq. 'THER_') then
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = .false.
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = .true.
        end if
        call ntfcma(zk16(icomp), zi(imate), aniso, ifon)
        if (aniso) then
            call jevech('PCAMASS', 'L', icamas)
        end if
    end if
!
    valQP = 0.0
    do kp = 1, FEQuadCell%nbQuadPoints
        tpg = FEEvalFuncScal(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
        dtpg = FEEvalGradVec(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
!
        if (zk16(icomp) (1:5) .eq. 'THER_') then
            call ntcomp(icomp, icamas, FECell%ndim, tpg, dtpg, &
                        FEQuadCell%points(1:3, kp), aniso, ifon, fluglo)
            flux(FECell%ndim*(kp-1)+1:FECell%ndim*(kp-1)+FECell%ndim) = -fluglo(1:FECell%ndim)
            valQP(1:3, kp) = fluglo
        else if (zk16(icomp) (1:5) .eq. 'SECH_') then
            tpsec = FEEvalFuncScal(FEBasis, sechf, FEQuadCell%points_param(1:3, kp))
            call rcdiff(zi(imate), zk16(icomp), tpsec, tpg, diff)
            valQP(1:3, kp) = diff*dtpg
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call FEStiffVecScal(FEQuadCell, FEBasis, valQP, resi)
!
    call writeVector("PRESIDU", nbDof, resi)
end subroutine
