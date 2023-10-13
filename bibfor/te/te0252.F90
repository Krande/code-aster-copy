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
subroutine te0252(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
!
    implicit none
#include "asterfort/addVecLumped.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_THER_RESI'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    type(FE_Cell) :: FECell, subFECell(4)
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer :: icodre(1)
    character(len=16) :: phenom
    real(kind=8) :: valQP(MAX_QP), tpgi, r8bid, funcEF(MAX_BS)
    real(kind=8) :: resi(MAX_BS), resi_sub(MAX_BS)
    real(kind=8) :: chal(1)
    integer :: kp, i, imate, icomp, itempi
    integer :: ifon(6)
    integer :: connec(4, 27), ise, nbSubCell, nbDof
    aster_logical :: aniso
    real(kind=8), pointer :: hydrgp(:) => null()
!
!-----------------------------------------------------------------------
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    nbDof = FEBasis%size
!
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PMATERC', 'L', imate)
!
    if (zk16(icomp) (1:5) .eq. 'THER_') then
!
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = ASTER_FALSE
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = ASTER_TRUE
        end if
        call ntfcma(zk16(icomp), zi(imate), aniso, ifon)
        if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
            call jevech('PHYDRPR', 'L', vr=hydrgp)
            call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                        ' ', 'THER_HYDR', 0, ' ', [r8bid], &
                        1, 'CHALHYDR', chal, icodre(1), 1)
        end if
    end if
!
    resi = 0.d0
!
    call FECell%splitLumped(nbSubCell, subFECell, connec)
    do ise = 1, nbSubCell
        call FEQuadCell%initCell(subFECell(ise), "MASS")
        call FEBasis%initCell(subFECell(ise))
!
        valQP = 0.0
        do kp = 1, FEQuadCell%nbQuadPoints
            tpgi = 0.d0
            funcEF = FEBasis%func(FEQuadCell%points_param(1:3, kp))
            do i = 1, FEBasis%size
                tpgi = tpgi+zr(itempi-1+connec(ise, i))*funcEF(i)
            end do
!
            if (zk16(icomp) (1:5) .eq. 'THER_') then
                call rcfode(ifon(1), tpgi, valQP(kp), r8bid)
                if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
                    valQP(kp) = valQP(kp)-chal(1)*hydrgp(kp)
                end if
            else if (zk16(icomp) (1:5) .eq. 'SECH_') then
                valQP(kp) = tpgi
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
        call FeMakeRhsScal(FEQuadCell, FEBasis, valQP, resi_sub)
!
        call addVecLumped(resi, resi_sub, ise, FEBasis%size, connec)
    end do
!
    call writeVector("PRESIDU", nbDof, resi)
!
end subroutine
