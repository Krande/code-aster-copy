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

subroutine te0244(option, nomte)
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
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/rcdiff.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/rcfode.h"
#include "asterfort/ntcomp.h"
#include "asterfort/utmess.h"
#include "asterfort/writeVector.h"
#include "asterfort/Behaviour_type.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
!
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_EVOLNI'
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
    character(len=16) :: phenom, rela_name
    real(kind=8) :: valQPM(MAX_QP), BGSEval(3, MAX_BS)
    real(kind=8) :: valQPMP(MAX_QP)
    real(kind=8) :: resi_f(MAX_BS), resi_m(MAX_BS), resi(MAX_BS)
    real(kind=8) :: resi_mp(MAX_BS), resi_p(MAX_BS), dfluxglo(3)
    real(kind=8) ::  deltat, theta, chal(1), diff, Kglo(3, 3)
    real(kind=8) :: beta, dbeta, tpg, dtpg(3), flux(3), sechpg, dsechpg(3)
    integer(kind=8) :: kp, imate, icamas, ifon(6), itemps, iret
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: lhyd, aniso
    real(kind=8), pointer :: tempi(:) => null()
    real(kind=8), pointer :: hydrpg(:) => null()
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    call FEQuadMass%initCell(FECell, "MASS")
    call FEQuadRigi%initCell(FECell, "RIGI")
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPER', 'L', vr=tempi)
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_name = compor(RELA_NAME)
!
    deltat = zr(itemps+1)
    theta = zr(itemps+2)

    if (rela_name(1:5) .eq. 'THER_') then
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = ASTER_FALSE
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = ASTER_TRUE
        end if
        call ntfcma(rela_name, zi(imate), aniso, ifon)
!       pour stopper le calcul si PCAMASS n'est pas disponible
        if (aniso) then
            call jevech('PCAMASS', 'L', icamas)
        end if
    end if
!
    resi_f = 0.d0
    do kp = 1, FEQuadRigi%nbQuadPoints
        BGSEval = FEBasis%grad(FEQuadRigi%points_param(1:3, kp), FEQuadRigi%jacob(1:3, 1:3, kp))
!
        if (rela_name(1:5) .eq. 'THER_') then
            tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuadRigi%points_param(1:3, kp))
            dtpg = FEEvalGradVec(FEBasis, tempi, FEQuadRigi%points_param(1:3, kp), BGSEval)
            call ntcomp(rela_name, FECell%ndim, tpg, dtpg, &
                        FEQuadRigi%points(1:3, kp), aniso, ifon, flux, Kglo, dfluxglo)
        else if (rela_name(1:5) .eq. 'SECH_') then
            sechpg = FEEvalFuncRScal(FEBasis, tempi, FEQuadRigi%points_param(1:3, kp))
            dsechpg = FEEvalGradVec(FEBasis, tempi, FEQuadRigi%points_param(1:3, kp), BGSEval)
            call rcvarc(' ', 'TEMP', '-', 'RIGI', kp, 1, tpg, iret)
            if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
            call rcdiff(zi(imate), rela_name, tpg, sechpg, diff)
            flux = diff*dsechpg
        else
            ASSERT(ASTER_FALSE)
        end if
        call FEStiffResiScalAdd(FEBasis, BGSEval, FEQuadRigi%weights(kp), flux, resi_f)
    end do
!
    if (rela_name(1:9) .eq. 'THER_HYDR') then
        lhyd = ASTER_TRUE
        call jevech('PHYDRPM', 'L', vr=hydrpg)
!
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', 'THER_HYDR', 0, ' ', [0.d0], &
                    1, 'CHALHYDR', chal, icodre, 1)
    else
        lhyd = ASTER_FALSE
    end if
!
    do kp = 1, FEQuadMass%nbQuadPoints
        tpg = FEEvalFuncRScal(FEBasis, tempi, FEQuadMass%points_param(1:3, kp))
        if (rela_name(1:5) .eq. 'THER_') then
            call rcfode(ifon(1), tpg, beta, dbeta)
            if (lhyd) then
                valQPMP(kp) = (dbeta*tpg-chal(1)*hydrpg(kp))
                valQPM(kp) = (beta-chal(1)*hydrpg(kp))
            else
                valQPMP(kp) = dbeta*tpg
                valQPM(kp) = beta
            end if
        else if (rela_name(1:5) .eq. 'SECH_') then
            valQPM(kp) = tpg
            valQPMP(kp) = tpg
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPM, resi_m)
    call FeMakeRhsScal(FEQuadMass, FEBasis, ValQPMP, resi_mp)
!
    resi = (theta-1.0d0)*resi_f+resi_m/deltat
    resi_p = (theta-1.0d0)*resi_f+resi_mp/deltat
!
    call writeVector('PVECTTR', FEBasis%size, resi)
    call writeVector('PVECTTI', FEBasis%size, resi_p)
!
end subroutine
