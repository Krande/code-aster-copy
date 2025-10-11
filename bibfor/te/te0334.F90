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
subroutine te0334(option, nomte)
!
    use FE_basis_module
    use FE_eval_module
    use FE_quadrature_module
    use FE_topo_module
    use BehaviourStrain_type
    use BehaviourStrain_module
!
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calcgr.h"
#include "asterfort/dmatmc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/granvi.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D
! Option: EPSP_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8), parameter :: ksp = 1, mxcmel = 54
    integer(kind=8), parameter :: nbsgm = 4
    aster_logical :: l_modi_cp
    real(kind=8) :: epsiPlas(mxcmel), epsiCreep(nbsgm)
    real(kind=8) :: epsiTota(6), epsiVarc(6), epsiMeca(6), sigma(4)
    integer(kind=8) :: nbVari, variIndxTemp, nbVariGranger
    real(kind=8) :: e, nu, c1, c2, trsig
    aster_logical :: l_creep, lTempInVari, lCplan, lDplan, lMetaLemaAni, lTHM
    integer(kind=8) :: jvMater, jvDisp, jvCompor, jvVari, jvSigm, jvEpsi, jvTime
    integer(kind=8) :: kpg, npg, ndim, nno, iSig, iEps, nbSig, nbEps
    integer(kind=8) :: jtab(7), iret
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
    real(kind=8) :: BGSEval(3, MAX_BS)
    real(kind=8) :: anglNaut(3), coorpg(3)
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword, relaComp, kit_comp_1, kit_comp_2
    type(All_Varc_Strain) :: allVarcStrain
    real(kind=8) :: tempkpg, time
    real(kind=8) :: d(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
    lCplan = lteatt('C_PLAN', 'OUI')
    lDplan = lteatt('D_PLAN', 'OUI')
    lTHM = lteatt('TYPMOD2', 'THM')
    l_modi_cp = .true.

! - Initialize a FE cell
    call FECell%init()
    ndim = FECell%ndim
    nno = FECell%nbnodes
    ASSERT(nno .le. 9)

! - Initialization of quadrature
    call FEQuad%initCell(FECell, fami)
    npg = FEQuad%nbQuadPoints
    ASSERT(npg .le. 9)

! - Initialization of basis functions
    call FEBasis%initCell(FECell)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Internal variables
    call jevech('PVARIGR', 'L', jvVari)
    call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, itab=jtab)
    nbVari = max(jtab(6), 1)*jtab(7)

! - Stresses
    call jevech('PCONTRR', 'L', jvSigm)
    nbSig = nbsigm()
    nbEps = nbSig
    ASSERT(nbSig .eq. nbsgm)

! - Orthotropic parameters: cannot use => zero
    anglNaut = 0.d0

! - Current time
    call jevech('PINSTR', 'L', jvTime)
    allVarcStrain%time = zr(jvTime)
    allVarcStrain%hasTime = ASTER_TRUE

! - Material parameters
    call jevech('PMATERC', 'L', jvMater)
    call get_elas_id(zi(jvMater), elasID, elasKeyword)
    if (elasID .ne. ELAS_ISOT) then
        call utmess('F', 'ELEMENTS6_2')
    end if

! - Behaviour
    call jevech('PCOMPOR', 'L', jvCompor)
    relaComp = zk16(jvCompor-1+RELA_NAME)
    kit_comp_1 = zk16(jvCompor-1+CREEP_NAME)
    kit_comp_2 = zk16(jvCompor-1+PLAS_NAME)
    lMetaLemaAni = (relaComp == "META_LEMA_ANI")

! - Stress plane warning
    if (lCplan) then
        if (relaComp .ne. 'VMIS_ISOT_LINE' .and. &
            relaComp(1:4) .ne. 'ELAS' .and. &
            relaComp .ne. 'VMIS_ISOT_TRAC') then
            call utmess('A', 'ELEMENTS6_3', sk=relaComp)
        end if
    end if

! - Detect Granger law (creep)
    if (relaComp(1:13) .ne. 'BETON_GRANGER' .and. &
        (relaComp .ne. 'KIT_DDI' .or. kit_comp_1(1:13) .ne. 'BETON_GRANGER')) then
        l_creep = ASTER_FALSE
    else
        call granvi("3D", nvi_=nbVariGranger)
        l_creep = ASTER_TRUE
    end if

! - Has temperature in internal state variable (== maximum)
    lTempInVari = ASTER_FALSE
    if (relaComp .eq. 'BETON_DOUBLE_DP') then
        variIndxTemp = 3
        lTempInVari = ASTER_TRUE
    else if (relaComp .eq. 'KIT_DDI') then
        if (kit_comp_2 .eq. 'BETON_DOUBLE_DP') then
            if (kit_comp_1(1:13) .eq. 'BETON_GRANGER') then
                variIndxTemp = nbVariGranger+3
                lTempInVari = ASTER_TRUE
            else
                call utmess('F', 'COMPOR5_76')
            end if
        end if
    end if

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Current coordinates of Gauss point
        coorpg = FEQuad%points_param(1:3, kpg)

! ----- Compute the gradient of the scalar basis
        BGSEval = FEBasis%grad(coorpg, FEQuad%jacob(1:3, 1:3, kpg))

! ----- Kinematic - Total strains
        epsiTota = FEEvalGradSymMat(FEBasis, zr(jvDisp), coorpg, BGSEval)
        epsiTota(4) = epsiTota(4)/sqrt(2.d0)

! ----- Detect external state variable
        call strainDetectVarc('+', lTHM, fami, kpg, ksp, &
                              allVarcStrain)

! ----- Get current temperature
        tempkpg = r8nnem()
        if (allVarcStrain%list(VARC_STRAIN_TEMP)%exist) then
            tempkpg = allVarcStrain%list(VARC_STRAIN_TEMP)%varcCurr(1)
        end if

! ----- Change temperature from internal variable (maximum) for BETON_DOUBLE_DP/BETON_GRANGER
        if (lTempInVari) then
            if (tempkpg .lt. zr(jvVari+(kpg-1)*nbVari+variIndxTemp-1)) then
                tempkpg = zr(jvVari+(kpg-1)*nbVari+variIndxTemp-1)
            end if
        end if

! ----- Set temperature
        allVarcStrain%hasTemp = ASTER_TRUE
        allVarcStrain%temp = tempkpg

! ----- Get elastic parameters (only isotropic elasticity)
        call get_elas_para(fami, zi(jvMater), '+', kpg, ksp, &
                           elasID, elasKeyword, &
                           time=time, temp=tempkpg, e_=e, nu_=nu)
        ASSERT(elasID .eq. ELAS_ISOT)

! ----- Compute non-mechanical strains (epsiVarc) for some external state variables
        call compVarcStrain(fami, '+', kpg, ksp, &
                            zi(jvMater), lMetaLemaAni, &
                            elasID, elasKeyword, &
                            allVarcStrain)
        call getVarcStrain('+', allVarcStrain, 6, epsiVarc)
        epsiVarc(4) = epsiVarc(4)/sqrt(2.d0)

! ----- Compute mechanical strains epsiMeca = epsiTota - epsiVarc
        epsiMeca = 0.d0
        epsiMeca(1:4) = epsiTota(1:4)-epsiVarc(1:4)

        if (lCplan) then
            call dmatmc(fami, zi(jvMater), time, '+', kpg, ksp, &
                        anglNaut, nbSig, d, l_modi_cp)
            epsiMeca(3) = -1.d0/d(3, 3)* &
                          (d(3, 1)*(epsiMeca(1)-epsiVarc(1))+ &
                           d(3, 2)*(epsiMeca(2)-epsiVarc(2))+ &
                           d(3, 4)*(epsiMeca(4)-2.d0*epsiVarc(4)))+epsiVarc(3)
        end if
        if (lDplan) then
            epsiMeca(3) = 0.d0
        end if

! ----- Compute creep strains
        epsiCreep = 0.d0
        if (l_creep) then
            call calcgr(kpg, nbSig, nbVari, zr(jvVari), nu, epsiCreep)
        end if

! ----- Compute stresses
        do iSig = 1, nbSig
            sigma(iSig) = zr(jvSigm+(kpg-1)*nbSig+iSig-1)
        end do
        if (lCplan) then
            trsig = sigma(1)+sigma(2)
        else
            trsig = sigma(1)+sigma(2)+sigma(3)
        end if

! ----- Compute plastic strains epsiPlas = epsi_tota - epsi_elas - epsiCreep
        c1 = (1.d0+nu)/e
        c2 = nu/e
        epsiPlas(nbEps*(kpg-1)+1) = epsiMeca(1)-(c1*sigma(1)-c2*trsig)-epsiCreep(1)
        epsiPlas(nbEps*(kpg-1)+2) = epsiMeca(2)-(c1*sigma(2)-c2*trsig)-epsiCreep(2)
        if (lCplan) then
            epsiPlas(nbEps*(kpg-1)+3) = -(epsiPlas(nbEps*(kpg-1)+1)+epsiPlas(nbEps*(kpg-1)+2))
        else
            epsiPlas(nbEps*(kpg-1)+3) = epsiMeca(3)-(c1*sigma(3)-c2*trsig)- &
                                        epsiCreep(3)
        end if
        epsiPlas(nbEps*(kpg-1)+4) = epsiMeca(4)-c1*sigma(4)-epsiCreep(4)
    end do

! - Plastic strain output
    call jevech('PDEFOPG', 'E', jvEpsi)
    do kpg = 1, npg
        do iEps = 1, nbEps
            zr(jvEpsi+nbEps*(kpg-1)+iEps-1) = epsiPlas(nbEps*(kpg-1)+iEps)
        end do
    end do
!
end subroutine
