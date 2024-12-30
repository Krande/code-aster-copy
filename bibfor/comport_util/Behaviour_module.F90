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
! aslint: disable=W1306
!
! ==================================================================================================
!
! Module for the management of integration of behaviour
!
! ==================================================================================================
!
module Behaviour_module
! ==================================================================================================
    use Behaviour_type
    use BehaviourMGIS_type
    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: behaviourOption, behaviourInit, behaviourInitPoint
    private :: initPara, initESVA
    public :: behaviourSetParaCell, behaviourSetParaPoin, setFromCompor
    public :: behaviourPrepESVAPoin, behaviourPrepStrain, behaviourPrepESVAExte
    public :: behaviourPrepESVAGeom, behaviourPrepModel, behaviourPredictionStress
    public :: behaviourCoorGauss
    private :: computeStrainESVA, computeStrainMeca
    private :: varcIsGEOM, isSolverIsExte
    private :: prepEltSize1, prepGradVelo
    private :: prepHygr, prepPtot, prepTemp, prepSech, prepHydr, prepEpsa, prepFields
    public :: getAsterVariableName, getMFrontVariableName
    private :: setFromOption, setFromCarcri
! ==================================================================================================
    private
#include "asterc/indik8.h"
#include "asterc/mgis_get_esvs.h"
#include "asterc/mgis_get_number_of_esvs.h"
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/isdeco.h"
#include "asterfort/jevech.h"
#include "asterfort/leverettIsotMeca.h"
#include "asterfort/matinv.h"
#include "asterfort/nmgeom.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! behaviourInit
!
! Initialisation of behaviour datastructure
!
! Out BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourInit(BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(out) :: BEHinteg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Initialization of datastructures'
        end if

! ----- Initialization of parameters
        call initPara(BEHinteg%behavPara)

! ----- Initialization of parameters for external state variables
        call initESVA(BEHinteg%behavESVA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initParaCell
!
! Initialisation of parameters
!
! Out BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine initPara(behavPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Para), intent(out) :: behavPara
!   ------------------------------------------------------------------------------------------------
!
        behavPara%timeCurr = r8nnem()
        behavPara%timePrev = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initESVA
!
! Initialization of parameters for external state variables
!
! Out behavESVA        : parameters for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine initESVA(behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(BehaviourESVA), intent(out) :: behavESVA
! ----- Locam
        aster_logical :: lGeomInESVA
        integer :: iField
!   ------------------------------------------------------------------------------------------------
!

! ----- Detect if geometry is in external state variables
        call varcIsGEOM(lGeomInESVA)
        behavESVA%lGeomInESVA = lGeomInESVA

! ----- For fields (default: None)
        do iField = 1, ESVA_FIELD_NBMAXI
            behavESVA%behavESVAField(iField)%exist = ASTER_FALSE
            behavESVA%behavESVAField(iField)%nbComp = 0
            behavESVA%behavESVAField(iField)%typeForStrain = ESVA_FIELD_TYPE_UNKW
            behavESVA%behavESVAField(iField)%valePrev = r8nnem()
            behavESVA%behavESVAField(iField)%valeCurr = r8nnem()
            behavESVA%behavESVAField(iField)%valeIncr = r8nnem()
            behavESVA%behavESVAField(iField)%valeScalRefe = r8nnem()
            behavESVA%behavESVAField(iField)%valeScalPrev = r8nnem()
            behavESVA%behavESVAField(iField)%valeScalIncr = r8nnem()
        end do

! ----- For geometric properties
        behavESVA%behavESVAGeom%elemSize1 = r8nnem()
        behavESVA%behavESVAGeom%gradVelo(9) = r8nnem()
        behavESVA%behavESVAGeom%coorElga(27, 3) = r8nnem()

! ----- For other properties
        behavESVA%behavESVAOther%rotpg(3*3) = r8nnem()
        behavESVA%behavESVAOther%tenscab = r8nnem()
        behavESVA%behavESVAOther%curvcab = r8nnem()
        behavESVA%behavESVAOther%nonloc(2) = r8nnem()
        behavESVA%behavESVAOther%r = r8nnem()
        behavESVA%behavESVAOther%time = r8nnem()
        behavESVA%behavESVAOther%hygrPrev = r8nnem()
        behavESVA%behavESVAOther%hygrIncr = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourSetParaCell
!
! Set parameters on a cell
!
! In  ldcDime          : dimension of physic for behaviour
! In  typmod           : finite element model
! In  option           : option to compute
! In  compor           : map for behaviour
! In  carcri           : parameters for comportment
! In  timePrev         : time at beginning of time step
! In  timeCurr         : time at end of time step
! In  fami             : Gauss family for integration point rule
! In  jvMaterCode      : adress for material parameters
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourSetParaCell(ldcDime, typmod, option, &
                                    compor, carcri, &
                                    timePrev, timeCurr, &
                                    fami, jvMaterCode, &
                                    BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: ldcDime
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: timePrev, timeCurr
        character(len=4), intent(in) :: fami
        integer, intent(in) :: jvMaterCode
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Local
        integer :: elasType, icodre
        character(len=16) :: elasKeyword
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Paramètres constants sur la cellule'
        end if

! ----- Set general parameters on cell
        BEHInteg%behavPara%ldcDime = ldcDime
        BEHInteg%behavPara%timePrev = timePrev
        BEHInteg%behavPara%timeCurr = timeCurr
        BEHInteg%behavPara%fami = fami

! ----- Set parameters for material properties
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Récupération de l élasticité'
        end if
        BEHInteg%behavPara%jvMaterCode = jvMaterCode
        call rccoma(jvMaterCode, 'ELAS', 0, elasKeyword, icodre)
        if (icodre .eq. 0) then
            call get_elas_id(jvMaterCode, elasType, elasKeyword)
            BEHInteg%behavPara%lElasIsMeta = (elasKeyword == 'ELAS_META')
            BEHInteg%behavPara%elasType = elasType
        end if
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Type d élasticité      : ', BEHInteg%behavPara%elasType
            WRITE (6, *) '<DEBUG>  Présence de métallurgie: ', BEHInteg%behavPara%lElasIsMeta
        end if

! ----- Set parameters from option
        call setFromOption(option, BEHinteg)

! ----- Set parameters from COMPOR map
        call setFromCompor(compor, BEHinteg)

! ----- Set parameters from CARCRI map
        if (BEHinteg%behavPara%lNonLinear) then
            call setFromCarcri(carcri, BEHinteg)
        end if

! ----- Set finite element model
        call behaviourPrepModel(typmod, BEHinteg)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromOption
!
! In  option           : option to compute
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromOption(option, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        type(Behaviour_Integ), intent(inout) :: BEHinteg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From OPTION'
        end if

        BEHinteg%behavPara%lImplex = option .eq. "RIGI_MECA_IMPLEX" .or. &
                                     option .eq. "RAPH_MECA_IMPLEX"
        BEHinteg%behavPara%lVari = L_VARI(option)
        BEHinteg%behavPara%lSigm = L_SIGM(option)
        BEHinteg%behavPara%lMatr = L_MATR(option)
        BEHinteg%behavPara%lPred = L_PRED(option)
        BEHinteg%behavPara%lNonLinear = option .ne. "FORC_NODA"
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromCompor
!
! In  compor           : map for behaviour
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromCompor(compor, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Locals
        character(len=16) :: defo_ldc, defo_comp, regu_visc, postIncr, mgisAddr
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From COMPOR'
        end if
        read (compor(DEFO_LDC), '(A16)') defo_ldc
        read (compor(DEFO), '(A16)') defo_comp
        read (compor(REGUVISC), '(A16)') regu_visc
        read (compor(POSTINCR), '(A16)') postIncr
        read (compor(MGIS_ADDR), '(A16)') mgisAddr
        BEHinteg%behavPara%lFiniteStrain = defo_comp .eq. 'SIMO_MIEHE' .or. &
                                           defo_comp .eq. 'GROT_GDEP'
        BEHinteg%behavPara%lGdefLog = defo_comp .eq. 'GDEF_LOG'
        BEHinteg%behavPara%lAnnealing = postIncr .eq. "REST_ECRO"
        BEHinteg%behavPara%lStrainMeca = defo_ldc .eq. 'MECANIQUE'
        BEHinteg%behavPara%lStrainAll = defo_ldc .eq. 'TOTALE'
        BEHinteg%behavPara%lStrainOld = defo_ldc .eq. 'OLD'
        BEHinteg%behavPara%lReguVisc = regu_visc .eq. 'REGU_VISC_ELAS'
        BEHinteg%behavESVA%behavESVAExte%mgisAddr = mgisAddr
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From COMPOR - lFiniteStrain: ', BEHinteg%behavPara%lFiniteStrain
            WRITE (6, *) '<DEBUG>  From COMPOR - lGdefLog: ', BEHinteg%behavPara%lGdefLog
            WRITE (6, *) '<DEBUG>  From COMPOR - lAnnealing: ', BEHinteg%behavPara%lAnnealing
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainMeca: ', BEHinteg%behavPara%lStrainMeca
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainAll: ', BEHinteg%behavPara%lStrainAll
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainOld: ', BEHinteg%behavPara%lStrainOld
            WRITE (6, *) '<DEBUG>  From COMPOR - lReguVisc: ', BEHinteg%behavPara%lReguVisc
            WRITE (6, *) '<DEBUG>  From COMPOR - mgisAddr: ', mgisAddr
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromCarcri
!
! In  carcri           : parameters for comportment
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromCarcri(carcri, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Locals
        aster_logical :: lMGIS, lUMAT
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From CARCRI'
        end if

! ----- External solvers ?
        call isSolverIsExte(carcri, lMGIS, lUMAT)
        BEHinteg%behavPara%lMGIS = lMGIS
        BEHinteg%behavPara%lUMAT = lUMAT
        BEHinteg%behavPara%lExteSolver = lUMAT .or. lMGIS

! ----- DEBUG
        if (LDC_PREP_DEBUG .eq. 1) then
        if (lMGIS) then
            WRITE (6, *) '<DEBUG>  External solver: MFront'
        elseif (lUMAT) then
            WRITE (6, *) '<DEBUG>  External solver: UMAT'
        else
            WRITE (6, *) '<DEBUG>  Internal solver'
        end if
        end if

! ----- Get list of external state variables from user
        call getListUserESVA(carcri, BEHinteg%behavESVA%tabcod)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepModel
!
! In  typmod           : finite element model
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepModel(typmod, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: typmod(2)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Locals
        integer :: lawOffset
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Paramètres du modèle'
        end if

        BEHinteg%behavPara%lStandardFE = typmod(2) .eq. ' '
        BEHinteg%behavPara%lTHM = typmod(2) .eq. 'THM'
        BEHinteg%behavPara%lCZM = typmod(2) .eq. 'ELEMJOIN'
        BEHinteg%behavPara%lGradVari = typmod(2) .eq. 'GRADVARI'
        BEHinteg%behavPara%lAxis = typmod(1) .eq. 'AXIS'
        BEHinteg%behavPara%lThreeDim = typmod(1) (1:2) .eq. '3D'
        BEHinteg%behavPara%lPlaneStrain = typmod(1) (1:6) .eq. 'D_PLAN'
        BEHinteg%behavPara%lPlaneStress = typmod(1) (1:6) .eq. 'C_PLAN'

        lawOffset = 0
        if (BEHinteg%behavPara%lImplex) then
            lawOffset = lawOffset+2000
        end if
        if (typmod(2) .eq. 'GDVARINO') then
            lawOffset = lawOffset+3000
        end if
        if (typmod(2) .eq. 'GRADSIGM') then
            lawOffset = lawOffset+4000
        end if
        if (typmod(2) .eq. 'GRADVARI') then
            lawOffset = lawOffset+6000
        end if
        if (typmod(2) .eq. 'EJ_HYME' .or. &
            typmod(2) .eq. 'ELEMJOIN' .or. typmod(2) .eq. 'INTERFAC') then
            lawOffset = lawOffset+7000
        end if
        BEHinteg%behavPara%lawIndexOffset = lawOffset
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepStrain
!
! Prepare input strains for the behaviour law
!    -> If defo_ldc = 'MECANIQUE', prepare mechanical strain
!    -> If defo_ldc = 'TOTALE' or 'OLD', keep total strain
!
! In  neps             : number of components of strains
! IO  epsm             : In : total strains at beginning of current step time
!                        Out : mechanical strains at beginning of current step time
! IO  deps             : In : increment of total strains during current step time
!                        Out : increment of mechanical strains during current step time
! In  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepStrain(neps, epsm, deps, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: neps
        real(kind=8), intent(inout) :: epsm(neps), deps(neps)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Local
        aster_logical :: lFiniteStrain, lPtot, lStrainMeca, lStrainAll
!   ------------------------------------------------------------------------------------------------
!
        if (ca_nbcvrc_ .ne. 0) then
            lStrainMeca = BEHinteg%behavPara%lStrainMeca
            lStrainAll = BEHinteg%behavPara%lStrainAll
            lFiniteStrain = BEHinteg%behavPara%lFiniteStrain
            lPtot = BEHinteg%behavESVA%behavESVAField(ESVA_FIELD_PTOT)%exist
            if (lStrainMeca .or. lPtot) then
                if (LDC_PREP_DEBUG .eq. 1) then
                    WRITE (6, *) '<DEBUG>  Présence de VARC avec nouveau système ou PTOT'
                end if
                ASSERT(.not. lFiniteStrain)
! ------------- Compute non-mechanic strains for some external state variables
                call computeStrainESVA(BEHinteg%behavESVA, &
                                       BEHinteg%behavPara%ldcDime, neps)
! ------------- Subtract to get mechanical strain epsm and deps become mechanical strains
                call computeStrainMeca(BEHinteg, neps, epsm, deps)
            end if
        else
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Absence de VARC'
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourSetParaPoin
!
! In  kpg              : index of quadrature point
! In  ksp              : index of "sub"-point (plates, pipes, beams, etc.)
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourSetParaPoin(kpg, ksp, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: kpg, ksp
        type(Behaviour_Integ), intent(inout) :: BEHinteg
!   ------------------------------------------------------------------------------------------------
!
        BEHinteg%behavPara%kpg = kpg
        BEHinteg%behavPara%ksp = ksp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourInitPoint
!
! Initialisation of behaviour datastructure - Special for SIMU_POINT_MAT
!
! In  relaComp         : RELATION in COMPORTemENT
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourInitPoint(relaComp, BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: relaComp
        type(Behaviour_Integ), intent(inout) :: BEHinteg
!   ------------------------------------------------------------------------------------------------
!

! ----- Set "real" zero
        BEHinteg%behavESVA%behavESVAGeom%coorElga = 0.d0

! ----- Special for GRAD_VELO
        if (BEHinteg%behavESVA%tabcod(GRADVELO) .eq. 1) then
            call utmess('A', 'COMPOR2_39')
            BEHinteg%behavESVA%behavESVAGeom%gradVelo = 0.d0
        end if

! ----- Special for ELTSIZE1
        if (BEHinteg%behavESVA%tabcod(ELTSIZE1) .eq. 1) then
            if (relaComp .ne. 'BETON_DOUBLE_DP') then
                call utmess('F', 'COMPOR2_12')
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListUserESVA
!
! Get list of external state variables from user (AFFE_VARC)
!
! In  carcri           : parameters for comportment
! Out tabcod           : list of integers to detect external state variables
!
!   ------------------------------------------------------------------------------------------------
    subroutine getListUserESVA(carcri, tabcod)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer, intent(out) :: tabcod(60)
! ----- Local
        integer :: jvariext1, jvariext2, variextecode(2)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Get list of external state variables from user'
        end if

        jvariext1 = nint(carcri(IVARIEXT1))
        jvariext2 = nint(carcri(IVARIEXT2))
        tabcod = 0
        variextecode(1) = jvariext1
        variextecode(2) = jvariext2
        call isdeco(variextecode, tabcod, 60)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAGeom
!
! Prepare external state variables - Geometry
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! IO  BEHinteg         : main object for managing the integration of behavior laws
! In  deplm            : displacements of nodes at beginning of time step
! In  ddepl            : displacements of nodes since beginning of time step
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAGeom(nno, npg, ndim, &
                                     jv_poids, jv_func, jv_dfunc, &
                                     geom, BEHinteg, &
                                     deplm_, ddepl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: nno, npg, ndim
        integer, intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        real(kind=8), optional, intent(in) :: deplm_(ndim, nno), ddepl_(ndim, nno)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for each element'
        end if

! ----- Compute element size 1
        if (BEHinteg%behavESVA%tabcod(ELTSIZE1) .eq. 1) then
            call prepEltSize1(nno, npg, ndim, &
                              jv_poids, jv_func, jv_dfunc, &
                              geom, BEHinteg%behavPara, &
                              BEHInteg%behavESVA%behavESVAGeom)
        end if

! ----- Compute gradient of velocity
        if (BEHinteg%behavESVA%tabcod(GRADVELO) .eq. 1) then
        if (.not. present(deplm_) .or. .not. present(ddepl_)) then
            call utmess('F', 'COMPOR2_26')
        end if
        call prepGradVelo(nno, npg, ndim, &
                          jv_poids, jv_func, jv_dfunc, &
                          geom, deplm_, ddepl_, &
                          BEHInteg%behavESVA%behavESVAGeom)
        end if

! ----- Coordinates of Gauss points (always)
        call behaviourCoorGauss(nno, npg, ndim, &
                                jv_func, geom, &
                                BEHInteg%behavESVA%behavESVAGeom)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAPoin
!
! Prepare external state variables at Gauss point
!
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAPoin(BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Local
        aster_logical :: lExteSolver, lStrainMeca, lhasInelasticStrains
        integer :: iField
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for current point'
        end if

! ----- Flags for external solvers
        lExteSolver = BEHinteg%behavPara%lExteSolver

! ----- Flag for strain decomposition
        lStrainMeca = BEHinteg%behavPara%lStrainMeca

! ----- Prepare hygrometry (not in AFFE_VARC)
        if (BEHinteg%behavESVA%tabcod(HYGR) .eq. 1) then
            call prepHygr(BEHinteg%behavPara%fami, &
                          BEHinteg%behavPara%kpg, &
                          BEHinteg%behavPara%ksp, &
                          BEHinteg%behavPara%jvMaterCode, &
                          BEHinteg%behavESVA)
        end if

! ----- Prepare external PTOT state variables (in AFFE_VARC)
        if (ca_nbcvrc_ .ne. 0) then
            call prepPtot(BEHinteg%behavPara%fami, &
                          BEHinteg%behavPara%kpg, &
                          BEHinteg%behavPara%ksp, &
                          BEHinteg%behavPara%jvMaterCode, &
                          BEHinteg%behavESVA)
        end if

! ----- Prepare external state variables from fields (in AFFE_VARC)
        if (ca_nbcvrc_ .ne. 0 .or. BEHinteg%behavPara%lTHM) then
            if (lExteSolver .or. lStrainMeca) then
                call prepFields(BEHinteg%behavPara%fami, &
                                BEHinteg%behavPara%kpg, &
                                BEHinteg%behavPara%ksp, &
                                BEHinteg%behavPara%jvMaterCode, &
                                BEHinteg%behavPara, &
                                BEHinteg%behavESVA)
            end if
        end if

! ----- Detect inelastic strains
        lhasInelasticStrains = ASTER_FALSE
        if (ca_nbcvrc_ .ne. 0) then
            do iField = 1, ESVA_FIELD_NBMAXI
                if (BEHinteg%behavESVA%behavESVAField(iField)%exist) then
                    lhasInelasticStrains = ASTER_TRUE
                    exit
                end if
            end do
        end if
        BEHinteg%behavESVA%lhasInelasticStrains = lhasInelasticStrains

! ----- Prepare other external state variables (For temperature: see preparation of fields)
        if (lExteSolver) then
            BEHinteg%behavESVA%behavESVAOther%time = BEHinteg%behavPara%timePrev
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMFrontVariableName
!
! Get MFront glossary or entry name from code_aster external state variable name
!
! In  exteNameAster    : code_aster external state variable name
!
! --------------------------------------------------------------------------------------------------
    character(len=64) function getMFrontVariableName(exteNameAster)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=64), intent(in) :: exteNameAster
! ----- Local
        integer :: iEsva
!   ------------------------------------------------------------------------------------------------
!
        getMFrontVariableName = exteNameAster
        do iEsva = 1, ESVA_EXTE_MGIS_NBMAXI
            if (fromAsterToMFront(1, iEsva) == exteNameAster) then
                getMFrontVariableName = fromAsterToMFront(2, iEsva)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! getAsterVariableName
!
! Get Aster name from MFront variable name
!
! In  exteNameMGIS     : MGIS external state variable name
!
! --------------------------------------------------------------------------------------------------
    character(len=8) function getAsterVariableName(exteNameMGIS)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=64), intent(in) :: exteNameMGIS
! ----- Local
        character(len=64) :: exteNameAster
        integer :: iEsva
!   ------------------------------------------------------------------------------------------------
!
        exteNameAster = exteNameMGIS
        do iEsva = 1, ESVA_EXTE_MGIS_NBMAXI
            if (fromAsterToMFront(2, iEsva) == exteNameMGIS) then
                exteNameAster = fromAsterToMFront(1, iEsva)
            end if
        end do
        getAsterVariableName = exteNameAster(1:8)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! logUndefinedVariable
!
! Log the lack of definition of external state variable
!
! In  exteNameAster    : variable name
!
! --------------------------------------------------------------------------------------------------
    subroutine logUndefinedVariable(exteNameAster)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: exteNameAster
! ----- Local
        character(len=64) :: exteNameAster64, MGISName
!   ------------------------------------------------------------------------------------------------
!
        exteNameAster64 = exteNameAster
        MGISName = getMFrontVariableName(exteNameAster64)
        if (MGISName .eq. exteNameAster64) then
            call utmess('F', 'COMPOR4_23', sk=exteNameAster)
        else
            call utmess('F', 'COMPOR4_75', nk=2, valk=[exteNameAster, MGISName])
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAExte
!
! Prepare external state variables for external solvers (UMAT/MFRONT)
!
! IO  BEHinteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAExte(BEHinteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHinteg
! ----- Local
        character(len=8), parameter :: nameUMAT(ESVA_EXTE_NBMAXI) = (/ &
                                       'SECH    ', 'HYDR    ', 'IRRA    ', &
                                       'NEUT1   ', 'NEUT2   ', 'CORR    ', &
                                       'ALPHPUR ', 'ALPHBETA'/)
        character(len=64) :: exteNameESVA(ESVA_EXTE_NBMAXI)
        character(len=64) :: exteNameMGIS
        character(len=8) :: exteNameAster
        real(kind=8) :: valePrev, valeCurr
        integer :: iESVA, nbESVA, indxField, iret
        aster_logical :: lMGIS, lUMAT, exist
        character(len=16) :: mgisAddr
        character(len=4) :: fami
        integer :: ksp, kpg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for external solvers'
        end if

! ----- Get parameters
        fami = BEHinteg%behavPara%fami
        kpg = BEHinteg%behavPara%kpg
        ksp = BEHinteg%behavPara%ksp
        lMGIS = BEHinteg%behavPara%lMGIS
        lUMAT = BEHinteg%behavPara%lUMAT

! ----- Get list of external state variables in external solvers
        nbESVA = 0
        if (lMGIS) then
            mgisAddr = BEHinteg%behavESVA%behavESVAExte%mgisAddr
            call mgis_get_number_of_esvs(mgisAddr, nbESVA)
            ASSERT(nbESVA .le. ESVA_EXTE_NBMAXI)
            call mgis_get_esvs(mgisAddr, exteNameESVA)
        elseif (lUMAT) then
            exteNameESVA = nameUMAT
            nbESVA = ESVA_EXTE_NBMAXI
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Default: all external states variables are scalar (not strain )
        BEHinteg%behavESVA%behavESVAExte%nbESVAScal = nbESVA

        if (LDC_PREP_DEBUG .eq. 1) then
            if (nbESVA .eq. 0) then
                WRITE (6, *) '<DEBUG> No external state variables defined in MFront/UMAT'
            else
                WRITE (6, *) '<DEBUG> Number of external state variables defined in MFront/UMAT:', &
                    nbESVA
            end if
        end if

! ----- Set values of ExternalStateVariables
        BEHinteg%behavESVA%behavESVAExte%scalESVAPrev = 0.d0
        BEHinteg%behavESVA%behavESVAExte%scalESVAIncr = 0.d0

        do iESVA = 1, nbESVA
! --------- Translate names of external state variable
            exteNameAster = getAsterVariableName(exteNameESVA(iESVA))
            exteNameMGIS = getMFrontVariableName(exteNameESVA(iESVA))

            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG> External state variable:', iESVA
                WRITE (6, *) '<DEBUG>  Aster name: ', exteNameAster
                if (lMGIS) then
                    WRITE (6, *) '<DEBUG>  MGIS name: ', exteNameMGIS
                end if
                if (lUMAT) then
                    WRITE (6, *) '<DEBUG>  UMAT name: ', exteNameESVA(iESVA)
                end if
            end if

            select case (exteNameMGIS)
            case ('GRADVELO')
                call utmess('F', 'COMPOR4_25', sk=exteNameAster)

            case ('ConcreteDrying')
                indxField = ESVA_FIELD_SECH
                exist = BEHInteg%behavESVA%behavESVAField(indxField)%exist
                if (exist) then
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalPrev
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalIncr
                else
                    if (.not. lUMAT) then
                        call logUndefinedVariable(exteNameAster)
                    end if
                end if

            case ('ConcreteHydration')
                indxField = ESVA_FIELD_HYDR
                exist = BEHInteg%behavESVA%behavESVAField(indxField)%exist
                if (exist) then
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalPrev
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalIncr
                else
                    if (BEHInteg%behavESVA%behavESVAOther%lHygr) then
                        BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = 0.d0
                        BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = 0.d0
                    else
                        if (.not. lUMAT) then
                            call logUndefinedVariable(exteNameAster)
                        end if
                    end if
                end if

            case ('Hygrometry')
                ASSERT(lMGIS)
                exist = BEHInteg%behavESVA%behavESVAOther%lHygr
                if (exist) then
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAOther%hygrPrev
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%behavESVA%behavESVAOther%hygrIncr
                else
                    call utmess('F', 'COMPOR4_26', sk=exteNameAster)
                end if

            CASE ('Temperature')
                ASSERT(lMGIS)
                indxField = ESVA_FIELD_TEMP
                BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                    BEHInteg%behavESVA%behavESVAField(indxField)%valeScalPrev
                BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                    BEHInteg%behavESVA%behavESVAField(indxField)%valeScalIncr

            case ('ElementSize')
                ASSERT(lMGIS)
                ASSERT(BEHInteg%behavESVA%behavESVAGeom%lElemSize1)
                BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                    BEHInteg%behavESVA%behavESVAGeom%elemSize1

            CASE ('ReferenceTemperature')
                ASSERT(lMGIS)
                indxField = ESVA_FIELD_TEMP
                exist = BEHInteg%behavESVA%behavESVAField(indxField)%exist
                if (exist) then
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalRefe
                else
                    call utmess('F', 'COMPOR4_26', sk=exteNameAster)
                end if

            CASE ('Time')
                ASSERT(lMGIS)
                BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                    BEHInteg%behavESVA%behavESVAOther%time

            CASE ('TEMP')
                ASSERT(lUMAT)
                indxField = ESVA_FIELD_TEMP
                exist = BEHInteg%behavESVA%behavESVAField(indxField)%exist
                if (exist) then
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalPrev
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%behavESVA%behavESVAField(indxField)%valeScalIncr
                else
                    call utmess('F', "COMPOR4_76")
                end if

            case default
                call rcvarc(' ', exteNameAster, '-', fami, kpg, ksp, valePrev, iret)
                if (iret .eq. 0) then
                    iret = 0
                    call rcvarc('F', exteNameAster, '+', fami, kpg, ksp, valeCurr, iret)
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = valePrev
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = valeCurr-valePrev
                else
                    if (.not. lUMAT) then
                        call logUndefinedVariable(exteNameAster)
                    end if
                end if
            end select
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Valeur précédente   : ', &
                    BEHinteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA)
                WRITE (6, *) '<DEBUG>  Valeur incrémentale : ', &
                    BEHinteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepPtot
!
! Prepare external state variable PTOT (specific)
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepPtot(fami, kpg, ksp, imate, &
                        behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: indxField = ESVA_FIELD_PTOT
        integer, parameter :: nbComp = 1
        integer :: iret
        aster_logical :: exist
        integer, parameter :: nbParaBiot = 1
        integer  :: codretBiot(nbParaBiot)
        real(kind=8) :: paraValeBiot(nbParaBiot)
        character(len=16), parameter :: paraNameBiot(nbParaBiot) = (/'BIOT_COEF'/)
        integer, parameter :: nbParaElas = 2
        integer  :: codretElas(nbParaElas)
        real(kind=8) :: paraValeElas(nbParaElas)
        character(len=16), parameter :: paraNameElas(nbParaElas) = (/'E ', 'NU'/)
        real(kind=8) :: ptotPrev, ptotCurr
        real(kind=8) :: ptotFieldPrev, ptotFieldCurr, ptotFieldIncr
        real(kind=8) :: biotp, biotm, em, ep, num, nup, troikm, troikp
!   ------------------------------------------------------------------------------------------------
!

! ----- Get values of external state variable
        exist = ASTER_FALSE
        ptotPrev = r8nnem()
        ptotCurr = r8nnem()
        iret = 0
        call rcvarc(' ', 'PTOT', '-', fami, kpg, ksp, &
                    ptotPrev, iret)
        if (iret .eq. 0) then
            iret = 0
            call rcvarc('F', 'PTOT', '+', fami, kpg, ksp, &
                        ptotCurr, iret)
            exist = ASTER_TRUE
        end if

! ----- Compute strains
        ptotFieldPrev = 0.d0
        ptotFieldCurr = 0.d0
        ptotFieldIncr = 0.d0
        if (exist) then
! --------- Get Biot coefficients
            call rcvalb(fami, kpg, ksp, &
                        '-', imate, ' ', 'THM_DIFFU', &
                        0, ' ', [0.d0], &
                        nbParaBiot, paraNameBiot, paraValeBiot, &
                        codretBiot, 1)
            if (codretBiot(1) .ne. 0) then
                paraValeBiot(1) = 0.d0
            end if
            biotm = paraValeBiot(1)
            call rcvalb(fami, kpg, ksp, &
                        '+', imate, ' ', 'THM_DIFFU', &
                        0, ' ', [0.d0], &
                        nbParaBiot, paraNameBiot, paraValeBiot, &
                        codretBiot, 1)
            if (codretBiot(1) .ne. 0) then
                paraValeBiot(1) = 0.d0
            end if
            biotp = paraValeBiot(1)

! --------- Get elastic coefficients
            call rcvalb(fami, kpg, ksp, &
                        '-', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbParaElas, paraNameElas, paraValeElas, &
                        codretElas, 1)
            if (codretElas(1) .ne. 0) then
                paraValeElas(1) = 0.d0
            end if
            if (codretElas(2) .ne. 0) then
                paraValeElas(2) = 0.d0
            end if
            em = paraValeElas(1)
            num = paraValeElas(2)
            call rcvalb(fami, kpg, ksp, &
                        '+', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbParaElas, paraNameElas, paraValeElas, &
                        codretElas, 1)
            if (codretElas(1) .ne. 0) then
                paraValeElas(1) = 0.d0
            end if
            if (codretElas(2) .ne. 0) then
                paraValeElas(2) = 0.d0
            end if
            ep = paraValeElas(1)
            nup = paraValeElas(2)
            troikp = ep/(1.d0-2.d0*nup)
            troikm = em/(1.d0-2.d0*num)
            ptotFieldPrev = (biotm/troikm)*ptotPrev
            ptotFieldCurr = (biotp/troikp)*ptotCurr
            ptotFieldIncr = ptotFieldCurr-ptotFieldPrev
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Prepare PTOT'
            end if
        end if

! ----- Save values
        behavESVA%behavESVAField(indxField)%exist = exist
        behavESVA%behavESVAField(indxField)%nbComp = nbComp
        behavESVA%behavESVAField(indxField)%typeForStrain = ESVA_FIELD_TYPE_VOLU
        behavESVA%behavESVAField(indxField)%valeCurr(1) = ptotFieldCurr
        behavESVA%behavESVAField(indxField)%valePrev(1) = ptotFieldPrev
        behavESVA%behavESVAField(indxField)%valeIncr(1) = ptotFieldIncr

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAField(indxField)%exist) then
                WRITE (6, *) '<DEBUG>  Values of PTOT: ', &
                    behavESVA%behavESVAField(indxField)%valePrev(1:nbComp), &
                    behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp)
            end if
        end if

!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepHygr
!
! Prepare external state variable HYGR (specific)
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepHygr(fami, kpg, ksp, imate, &
                        behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: nbPara = 1
        integer  :: codret(nbPara)
        real(kind=8) :: paraVale(nbPara)
        character(len=16), parameter :: paraName(nbPara) = (/'FONC_DESORP'/)
        character(len=16) :: phenom
        real(kind=8) :: funcDesorpPrev, funcDesorpCurr
        aster_logical :: exist
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare HYGR'
        end if

! ----- Get parameters
        funcDesorpPrev = 0.d0
        funcDesorpCurr = 0.d0
        exist = ASTER_FALSE
        call rccoma(imate, 'BETON_DESORP', 0, phenom, codret(1))
        if (codret(1) .ne. 0) then
            call utmess('F', 'COMPOR2_94')
        end if
        exist = ASTER_TRUE
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                    nbPara, paraName, paraVale, codret, 0)
        if (codret(1) .eq. 0) then
            funcDesorpPrev = paraVale(1)
            call rcvalb(fami, kpg, ksp, '+', imate, &
                        ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, codret, 0)
            ASSERT(codret(1) .eq. 0)
            funcDesorpCurr = paraVale(1)
        else
            !   leverett isotherm
            call leverettIsotMeca(fami, kpg, ksp, imate, funcDesorpPrev, &
                                  funcDesorpCurr)
        end if

! ----- Save values
        if (exist) then
            behavESVA%behavESVAOther%lHygr = exist
            behavESVA%behavESVAOther%hygrPrev = funcDesorpPrev
            behavESVA%behavESVAOther%hygrIncr = funcDesorpCurr-funcDesorpPrev
        end if

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAOther%lHygr) then
                WRITE (6, *) '<DEBUG>  Values of HYGR: ', &
                    behavESVA%behavESVAOther%hygrPrev, &
                    behavESVA%behavESVAOther%hygrIncr
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepTemp
!
! Prepare external state variable TEMP
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! In  behavPara        : parameters for integration of behaviour
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepTemp(fami, kpg, ksp, imate, &
                        behavPara, behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(Behaviour_Para), intent(in) :: behavPara
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: indxField = ESVA_FIELD_TEMP
        integer :: iret, nbComp
        aster_logical :: exist
        real(kind=8) :: tempRefe, tempPrev, tempCurr
        real(kind=8) :: tempFieldIncr(ESVA_FIELD_TEMP_NBCMP)
        real(kind=8) :: tempFieldPrev(ESVA_FIELD_TEMP_NBCMP)
        real(kind=8) :: tempFieldCurr(ESVA_FIELD_TEMP_NBCMP)
        real(kind=8) :: therIsotPrev, therAnisPrev(ESVA_FIELD_TEMP_NBCMP), therMetaPrev
        real(kind=8) :: therIsotCurr, therAnisCurr(ESVA_FIELD_TEMP_NBCMP), therMetaCurr
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare TEMP'
        end if

! ----- Get thermal values
        exist = ASTER_FALSE
        tempPrev = r8nnem()
        tempCurr = r8nnem()
        tempRefe = r8nnem()
        iret = 0
        call verift(fami, kpg, ksp, '-', imate, &
                    epsth_=therIsotPrev, &
                    epsth_anis_=therAnisPrev, &
                    epsth_meta_=therMetaPrev, &
                    temp_prev_=tempPrev, &
                    temp_refe_=tempRefe, &
                    iret_=iret)
        if (iret .ne. 0) then
            tempPrev = 0.d0
        end if
        iret = 0
        call verift(fami, kpg, ksp, '+', imate, &
                    epsth_=therIsotCurr, &
                    epsth_anis_=therAnisCurr, &
                    epsth_meta_=therMetaCurr, &
                    temp_curr_=tempCurr, &
                    iret_=iret)
        if (iret .eq. 0) then
            exist = ASTER_TRUE
        else
            tempCurr = 0.d0
        end if

! ----- Temperature as scalar (no strain)
        if (exist) then
            behavESVA%behavESVAField(indxField)%valeScalRefe = tempRefe
            behavESVA%behavESVAField(indxField)%valeScalPrev = tempPrev
            behavESVA%behavESVAField(indxField)%valeScalIncr = tempCurr-tempPrev
        end if

! ----- Number of components
        nbComp = 0
        if (exist) then
            if (behavPara%elasType .eq. ELAS_ISOT) then
                nbComp = 1
            elseif (behavPara%elasType .eq. ELAS_ISTR) then
                nbComp = 2
            elseif (behavPara%elasType .eq. ELAS_ORTH) then
                nbComp = 3
            end if
        end if
        ASSERT(nbComp .le. ESVA_FIELD_TEMP_NBCMP)

! ----- Compute THER strains
        tempFieldPrev = 0.d0
        tempFieldCurr = 0.d0
        tempFieldIncr = 0.d0
        if (exist) then
            if (behavPara%lElasIsMeta) then
                ASSERT(behavPara%elasType .eq. ELAS_ISOT)
                ASSERT(nbComp .eq. 1)
                tempFieldPrev(1:nbComp) = therMetaPrev
                tempFieldCurr(1:nbComp) = therMetaCurr
                tempFieldIncr(1:nbComp) = therMetaCurr-therMetaPrev
            else
                if (behavPara%elasType == ELAS_ISOT) then
                    ASSERT(nbComp .eq. 1)
                    tempFieldIncr(1:nbComp) = therIsotCurr-therIsotPrev
                    tempFieldPrev(1:nbComp) = therIsotPrev
                    tempFieldCurr(1:nbComp) = therIsotCurr
                elseif (behavPara%elasType == ELAS_ORTH .or. &
                        behavPara%elasType == ELAS_ISTR) then
                    tempFieldIncr(1:nbComp) = therAnisCurr(1:nbComp)-therAnisPrev(1:nbComp)
                    tempFieldPrev(1:nbComp) = therAnisPrev(1:nbComp)
                    tempFieldCurr(1:nbComp) = therAnisCurr(1:nbComp)
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        end if

! ----- Save values
        behavESVA%behavESVAField(indxField)%exist = exist
        behavESVA%behavESVAField(indxField)%nbComp = nbComp
        behavESVA%behavESVAField(indxField)%typeForStrain = ESVA_FIELD_TYPE_VOLU
        behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp) = tempFieldCurr(1:nbComp)
        behavESVA%behavESVAField(indxField)%valePrev(1:nbComp) = tempFieldPrev(1:nbComp)
        behavESVA%behavESVAField(indxField)%valeIncr(1:nbComp) = tempFieldIncr(1:nbComp)

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAField(indxField)%exist) then
                WRITE (6, *) '<DEBUG>  Values of TEMP (prev): ', &
                    behavESVA%behavESVAField(indxField)%valePrev(1:nbComp)
                WRITE (6, *) '<DEBUG>  Values of TEMP (curr): ', &
                    behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepSech
!
! Prepare external state variable SECH
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepSech(fami, kpg, ksp, imate, &
                        behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: indxField = ESVA_FIELD_SECH
        integer, parameter :: nbComp = 1
        integer :: iret
        aster_logical :: exist
        integer, parameter :: nbPara = 1
        real(kind=8) :: paraVale(nbPara)
        character(len=16), parameter :: paraName(nbPara) = ('K_DESSIC')
        integer :: codret(nbPara)
        real(kind=8) :: sechRefe, sechPrev, sechCurr
        real(kind=8) :: sechFieldIncr, sechFieldPrev, sechFieldCurr
        real(kind=8) :: kdessm, kdessp
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare SECH'
        end if

! ----- Get values of external state variable
        exist = ASTER_FALSE
        sechPrev = r8nnem()
        sechCurr = r8nnem()
        sechRefe = r8nnem()
        iret = 0
        call rcvarc(' ', 'SECH', '-', fami, kpg, ksp, sechPrev, iret)
        if (iret .eq. 0) then
            iret = 0
            call rcvarc('F', 'SECH', '+', fami, kpg, ksp, sechCurr, iret)
            exist = ASTER_TRUE
        end if
        iret = 0
        call rcvarc(' ', 'SECH', 'REF', fami, kpg, ksp, sechRefe, iret)
        if (iret .ne. 0) then
            sechRefe = 0.d0
        end if

! ----- Drying as scalar (no strain)
        if (exist) then
            behavESVA%behavESVAField(indxField)%valeScalRefe = sechRefe
            behavESVA%behavESVAField(indxField)%valeScalPrev = sechPrev
            behavESVA%behavESVAField(indxField)%valeScalIncr = sechCurr-sechPrev
        end if

! ----- Compute strains
        sechFieldPrev = 0.d0
        sechFieldCurr = 0.d0
        sechFieldIncr = 0.d0
        if (exist) then
            call rcvalb(fami, kpg, ksp, &
                        '-', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, &
                        codret, 1)
            kdessm = paraVale(1)
            call rcvalb(fami, kpg, ksp, &
                        '+', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, &
                        codret, 1)
            kdessp = paraVale(1)
            sechFieldPrev = -kdessm*(sechRefe-sechPrev)
            sechFieldCurr = -kdessp*(sechRefe-sechCurr)
            sechFieldIncr = sechFieldCurr-sechFieldPrev
        end if

! ----- Save values
        behavESVA%behavESVAField(indxField)%exist = exist
        behavESVA%behavESVAField(indxField)%nbComp = nbComp
        behavESVA%behavESVAField(indxField)%typeForStrain = ESVA_FIELD_TYPE_VOLU
        behavESVA%behavESVAField(indxField)%valeCurr(1) = sechFieldCurr
        behavESVA%behavESVAField(indxField)%valePrev(1) = sechFieldPrev
        behavESVA%behavESVAField(indxField)%valeIncr(1) = sechFieldIncr

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAField(indxField)%exist) then
                WRITE (6, *) '<DEBUG>  Values of SECH: ', &
                    behavESVA%behavESVAField(indxField)%valePrev(1:nbComp), &
                    behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepHydr
!
! Prepare external state variable HYDR
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepHydr(fami, kpg, ksp, imate, &
                        behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: indxField = ESVA_FIELD_HYDR
        integer, parameter :: nbComp = 1
        integer :: iret
        aster_logical :: exist
        integer, parameter :: nbPara = 1
        real(kind=8) :: paraVale(nbPara)
        character(len=16), parameter :: paraName(nbPara) = ('B_ENDOGE')
        integer :: codret(nbPara)
        real(kind=8) :: hydrPrev, hydrCurr
        real(kind=8) :: hydrFieldIncr, hydrFieldPrev, hydrFieldCurr
        real(kind=8) :: bendom, bendop
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare HYDR'
        end if

! ----- Get values of external state variables
        exist = ASTER_FALSE
        hydrPrev = r8nnem()
        hydrCurr = r8nnem()
        iret = 0
        call rcvarc(' ', 'HYDR', '-', fami, kpg, ksp, hydrPrev, iret)
        if (iret .eq. 0) then
            iret = 0
            call rcvarc('F', 'HYDR', '+', fami, kpg, ksp, hydrCurr, iret)
            exist = ASTER_TRUE
        end if

! ----- Hydratation as scalar (no strain)
        if (exist) then
            behavESVA%behavESVAField(indxField)%valeScalPrev = hydrPrev
            behavESVA%behavESVAField(indxField)%valeScalIncr = hydrCurr-hydrPrev
        end if

! ----- Compute strains
        hydrFieldPrev = 0.d0
        hydrFieldCurr = 0.d0
        hydrFieldIncr = 0.d0
        if (exist) then
            call rcvalb(fami, kpg, ksp, &
                        '-', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, &
                        codret, 1)
            bendom = paraVale(1)
            call rcvalb(fami, kpg, ksp, &
                        '+', imate, ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, &
                        codret, 1)
            bendop = paraVale(1)
            hydrFieldPrev = -bendom*hydrPrev
            hydrFieldCurr = -bendop*hydrCurr
            hydrFieldIncr = hydrFieldCurr-hydrFieldPrev
        end if

! ----- Save values
        behavESVA%behavESVAField(indxField)%exist = exist
        behavESVA%behavESVAField(indxField)%nbComp = nbComp
        behavESVA%behavESVAField(indxField)%typeForStrain = ESVA_FIELD_TYPE_VOLU
        behavESVA%behavESVAField(indxField)%valeCurr = hydrFieldCurr
        behavESVA%behavESVAField(indxField)%valePrev = hydrFieldPrev
        behavESVA%behavESVAField(indxField)%valeIncr = hydrFieldIncr

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAField(indxField)%exist) then
                WRITE (6, *) '<DEBUG>  Values of HYDR: ', &
                    behavESVA%behavESVAField(indxField)%valePrev(1:nbComp), &
                    behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepEpsa
!
! Prepare external state variable EPSA
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  behavPara        : parameters for integration of behaviour
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepEpsa(fami, kpg, ksp, &
                        behavPara, behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp
        type(Behaviour_Para), intent(in) :: behavPara
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer, parameter :: indxField = ESVA_FIELD_EPSA
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        integer :: iret, iComp, nbComp
        aster_logical :: exist
        character(len=6), parameter :: epsaName(ESVA_FIELD_EPSA_NBCMP) = &
                                       (/'EPSAXX', 'EPSAYY', 'EPSAZZ', &
                                         'EPSAXY', 'EPSAXZ', 'EPSAYZ'/)
        real(kind=8) :: epsaPrev(ESVA_FIELD_EPSA_NBCMP)
        real(kind=8) :: epsaCurr(ESVA_FIELD_EPSA_NBCMP)
        real(kind=8) :: epsaIncr(ESVA_FIELD_EPSA_NBCMP)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare EPSA'
        end if

! ----- Number of components
        if (behavPara%ldcDime .eq. 1) then
            nbComp = 1
        else
            nbComp = 6
        end if
        ASSERT(nbComp .le. ESVA_FIELD_EPSA_NBCMP)

! ----- Get values of external state variable
        exist = ASTER_FALSE
        epsaPrev = r8nnem()
        epsaCurr = r8nnem()
        iret = 0
        do iComp = 1, nbComp
            iret = 0
            call rcvarc(' ', epsaName(iComp), '-', fami, kpg, ksp, &
                        epsaPrev(iComp), iret)
            if (iret .ne. 0) then
                epsaPrev(iComp) = 0.d0
            else
                exist = ASTER_TRUE
            end if
            iret = 0
            call rcvarc(' ', epsaName(iComp), '+', fami, kpg, ksp, &
                        epsaCurr(iComp), iret)
            if (iret .ne. 0) then
                epsaCurr(iComp) = 0.d0
            else
                exist = ASTER_TRUE
            end if
            epsaIncr(iComp) = epsaCurr(iComp)-epsaPrev(iComp)
        end do

! ----- Nondiagonal terms of EPSA are rescaled with rac2
        do iComp = 4, nbComp
            epsaPrev(iComp) = epsaPrev(iComp)*rac2
            epsaCurr(iComp) = epsaCurr(iComp)*rac2
            epsaIncr(iComp) = epsaIncr(iComp)*rac2
        end do

! ----- Save values
        behavESVA%behavESVAField(indxField)%exist = exist
        behavESVA%behavESVAField(indxField)%nbComp = nbComp
        behavESVA%behavESVAField(indxField)%typeForStrain = ESVA_FIELD_TYPE_COMP
        behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp) = epsaCurr(1:nbComp)
        behavESVA%behavESVAField(indxField)%valePrev(1:nbComp) = epsaPrev(1:nbComp)
        behavESVA%behavESVAField(indxField)%valeIncr(1:nbComp) = epsaIncr(1:nbComp)

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAField(indxField)%exist) then
                WRITE (6, *) '<DEBUG>  Values of EPSA: ', &
                    behavESVA%behavESVAField(indxField)%valePrev(1:nbComp), &
                    behavESVA%behavESVAField(indxField)%valeCurr(1:nbComp)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepFields
!
! Prepare fields for external state variables (general)
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  imate            : coded material address
! In  behavPara        : parameters for integration of behaviour
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepFields(fami, kpg, ksp, imate, &
                          behavPara, behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        integer, intent(in) :: kpg, ksp, imate
        type(Behaviour_Para), intent(in) :: behavPara
        type(BehaviourESVA), intent(inout) :: behavESVA
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare fields for external state variables'
        end if

! ----- Prepare external state variable TEMP
        call prepTemp(fami, kpg, ksp, imate, &
                      behavPara, behavESVA)

! ----- Prepare external state variable SECH
        call prepSech(fami, kpg, ksp, imate, &
                      behavESVA)

! ----- Prepare external state variable HYDR
        call prepHydr(fami, kpg, ksp, imate, &
                      behavESVA)

! ----- Prepare external state variable EPSA
        call prepEpsa(fami, kpg, ksp, &
                      behavPara, behavESVA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPredictionStress
!
! Compute the contribution of the "thermal stress increment" in prediction
!
! in  behavESVA        : parameters for external state variables
! in  dsidep           : tangent operator in prediction
! IO  sig              : stress (Taylor term of order 0)
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPredictionStress(behavESVA, dsidep, sig)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(BehaviourESVA), intent(in) :: behavESVA
        real(kind=8), intent(in) :: dsidep(:, :)
        real(kind=8), intent(inout) :: sig(:)
! ----- Local
        integer :: ndimsi
!   ------------------------------------------------------------------------------------------------
!
        if (ca_nbcvrc_ .ne. 0) then
            ndimsi = size(sig)
            ASSERT(ndimsi .le. 6)
            ASSERT(size(dsidep, 1) .ge. ndimsi)
            ASSERT(size(dsidep, 2) .ge. ndimsi)
            sig = sig-matmul(dsidep(1:ndimsi, 1:ndimsi), behavESVA%depsi_varc(1:ndimsi))
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Prediction stress: ', sig
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepEltSize1
!
! Compute external state variables - Geometry: size of element (ELTSIZE1)
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! In  behavPara        : parameters for integration of behaviour
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine prepEltSize1(nno, npg, ndim, &
                            jv_poids, jv_func, jv_dfunc, &
                            geom, behavPara, &
                            behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: nno, npg, ndim
        integer, intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(Behaviour_Para), intent(in) :: behavPara
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer :: kpg, iComp, i
        real(kind=8) :: lc, dfdx(MT_NNOMAX3D), dfdy(MT_NNOMAX3D), dfdz(MT_NNOMAX3D), poids, r
        real(kind=8) :: volume, surfac
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nno .le. MT_NNOMAX3D)
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute ELTSIZE1'
        end if
!
        if (behavPara%lThreeDim) then
            volume = 0.d0
            do kpg = 1, npg
                call dfdm3d(nno, kpg, jv_poids, jv_dfunc, geom, &
                            poids, dfdx, dfdy, dfdz)
                volume = volume+poids
            end do
            if (npg .ge. 9) then
                lc = volume**0.33333333333333d0
            else
                lc = rac2*volume**0.33333333333333d0
            end if
        elseif (behavPara%lAxis .or. behavPara%lPlaneStrain) then
            surfac = 0.d0
            do kpg = 1, npg
                iComp = (kpg-1)*nno
                call dfdm2d(nno, kpg, jv_poids, jv_dfunc, geom, &
                            poids, dfdx, dfdy)
                if (behavPara%lAxis) then
                    r = 0.d0
                    do i = 1, nno
                        r = r+geom(1, i)*zr(jv_func+i+iComp-1)
                    end do
                    poids = poids*r
                end if
                surfac = surfac+poids
            end do
            if (npg .ge. 5) then
                lc = surfac**0.5d0
            else
                lc = rac2*surfac**0.5d0
            end if
        elseif (behavPara%lPlaneStress) then
            lc = r8nnem()
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save values
        behavESVAGeom%lElemSize1 = ASTER_TRUE
        behavESVAGeom%elemSize1 = lc
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepGradVelo
!
! Compute external state variables - Geometry: gradient of velocity (GRADVELO)
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! In  deplm            : displacements of nodes at beginning of time step
! In  ddepl            : displacements of nodes since beginning of time step
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine prepGradVelo(nno, npg, ndim, &
                            jv_poids, jv_func, jv_dfunc, &
                            geom, deplm, ddepl, &
                            behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: nno, npg, ndim
        integer, intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        real(kind=8), intent(in) :: deplm(ndim, nno), ddepl(ndim, nno)
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer :: nddl, kpg, i, j
        real(kind=8) :: l(3, 3), fmm(3, 3), df(3, 3), f(3, 3), r8bid, r
        real(kind=8) :: deplp(3, MT_NNOMAX3D), geomm(3, MT_NNOMAX3D), epsbid(6)
        real(kind=8) :: dfdi(nno, 3)
        real(kind=8), parameter :: id(9) = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        nddl = ndim*nno
        behavESVAGeom%gradVelo = 0.d0
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute GRADVELO'
        end if
!
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, geom, b_incx, geomm, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deplm, b_incx, geomm, &
                   b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deplm, b_incx, deplp, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, ddepl, b_incx, deplp, &
                   b_incy)
        do kpg = 1, npg
            call nmgeom(ndim, nno, .false._1, .true._1, geom, &
                        kpg, jv_poids, jv_func, jv_dfunc, deplp, &
                        .true._1, r8bid, dfdi, f, epsbid, &
                        r)
            call nmgeom(ndim, nno, .false._1, .true._1, geomm, &
                        kpg, jv_poids, jv_func, jv_dfunc, ddepl, &
                        .true._1, r8bid, dfdi, df, epsbid, &
                        r)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, id, b_incx, df, &
                       b_incy)
            call matinv('S', 3, f, fmm, r8bid)
            l = matmul(df, fmm)
            do i = 1, 3
                do j = 1, 3
                    behavESVAGeom%gradVelo(3*(i-1)+j) = l(i, j)
                end do
            end do
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourCoorGauss
!
! Compute external state variables - Geometry: coordinates of Gauss points
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_func          : JEVEUX adress for shape functions
! In  geom             : initial coordinates of nodes
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourCoorGauss(nno, npg, ndim, &
                                  jv_func, geom, &
                                  behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: nno, npg, ndim
        integer, intent(in) :: jv_func
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer :: i, iComp, kpg
!   ------------------------------------------------------------------------------------------------
!
        behavESVAGeom%coorElga = 0.d0
        ASSERT(npg .le. ESVA_GEOM_NBMAXI)
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute coordinates of Gauss points'
        end if
!
        do kpg = 1, npg
            do i = 1, ndim
                do iComp = 1, nno
                    behavESVAGeom%coorElga(kpg, i) = &
                        behavESVAGeom%coorElga(kpg, i)+ &
                        geom(i, iComp)*zr(jv_func-1+nno*(kpg-1)+iComp)
                end do
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! computeStrainESVA
!
! Pre-compute non-mechanical strains from external state variables
!
! IO  behavESVA        : parameters for External State Variables
! In  ldcDime          : number of components for volumic strains
! In  neps             : number of components of strains
!
! --------------------------------------------------------------------------------------------------
    subroutine computeStrainESVA(behavESVA, ldcDime, neps)
!   -----------------------------------------------------------------------------------------------
! ----- Parameters
        type(BehaviourESVA), intent(inout) :: behavESVA
        integer, intent(in) :: ldcDime
        integer, intent(in) :: neps
! ----- Local
        integer :: iComp, iField, iDime
        aster_logical :: exist, hasAnelastiStrains
        integer :: nbComp, typeForStrain
        real(kind=8) :: valePrev(ESVA_FIELD_NBCMPMAXI)
        real(kind=8) :: valeIncr(ESVA_FIELD_NBCMPMAXI)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Pre-compute non-mechanical strains from external state variables'
        end if

        behavESVA%depsi_varc = 0.d0
        behavESVA%epsi_varc = 0.d0
        hasAnelastiStrains = ASTER_FALSE

        do iField = 1, ESVA_FIELD_NBMAXI
! --------- Get current field
            exist = behavESVA%behavESVAField(iField)%exist
            nbComp = behavESVA%behavESVAField(iField)%nbComp
            valePrev = behavESVA%behavESVAField(iField)%valePrev
            valeIncr = behavESVA%behavESVAField(iField)%valeIncr
            typeForStrain = behavESVA%behavESVAField(iField)%typeForStrain

! --------- Add field if exist
            if (exist) then
                if (LDC_PREP_DEBUG .eq. 1) then
                    WRITE (6, *) '<DEBUG>  Field: ', iField, nbComp, &
                        valePrev(1:nbComp), valeIncr(1:nbComp)
                end if
                if (typeForStrain .eq. ESVA_FIELD_TYPE_COMP) then
                    if (nbComp .ne. neps) then
                        call utmess("F", "COMPOR4_77")
                    end if
                    if (LDC_PREP_DEBUG .eq. 1) then
                        WRITE (6, *) '<DEBUG> Déformation donnée par composante'
                    end if
                    do iComp = 1, nbComp
                        behavESVA%epsi_varc(iComp) = behavESVA%epsi_varc(iComp)+ &
                                                     valePrev(iComp)
                        behavESVA%depsi_varc(iComp) = behavESVA%depsi_varc(iComp)+ &
                                                      valeIncr(iComp)
                    end do
                else if (typeForStrain .eq. ESVA_FIELD_TYPE_VOLU) then
                    ASSERT(ldcDime .eq. 2 .or. ldcDime .eq. 3)
                    if (LDC_PREP_DEBUG .eq. 1) then
                        WRITE (6, *) '<DEBUG> Déformation volumique'
                    end if
                    if (nbComp .eq. 1) then
                        iComp = 1
                    else
                        call utmess("F", "COMPOR4_79")
                    end if
                    do iDime = 1, 3
                        behavESVA%epsi_varc(iDime) = behavESVA%epsi_varc(iDime)+ &
                                                     valePrev(iComp)
                        behavESVA%depsi_varc(iDime) = behavESVA%depsi_varc(iDime)+ &
                                                      valeIncr(iComp)
                    end do
                else
                    WRITE (6, *) "typeForStrain: ", typeForStrain
                    ASSERT(ASTER_FALSE)
                end if
            end if
        end do

! ----- DEBUG
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Strains from external state variables - Prev: ', &
                neps, behavESVA%epsi_varc(1:neps)
            WRITE (6, *) '<DEBUG>  Strains from external state variables - Incr: ', &
                neps, behavESVA%depsi_varc(1:neps)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! computeStrainMeca
!
! Prepare strains (substracting "thermic" strains to total strains to get mechanical part)
!
! In  BEHinteg         : main object for managing the integration of behavior laws
! In  neps             : number of components of strains
! IO  epsm             : In : total strains at beginning of current step time
!                        Out : mechanical strains at beginning of current step time
! IO  deps             : In : increment of total strains during current step time
!                        Out : increment of mechanical strains during current step time
!
! --------------------------------------------------------------------------------------------------
    subroutine computeStrainMeca(BEHinteg, neps, epsm, deps)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer, intent(in) :: neps
        real(kind=8), intent(inout) :: epsm(neps), deps(neps)
! ----- Local
        real(kind=8) :: stran(12), dstran(12)
        integer :: nepu
        aster_logical :: lGradVari, lCZM, lEpsa
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Calcul des déformations mécaniques'
        end if
        lCZM = BEHinteg%behavPara%lCZM
        lGradVari = BEHinteg%behavPara%lGradVari
        lEpsa = BEHinteg%behavESVA%behavESVAField(ESVA_FIELD_EPSA)%exist
        dstran = 0.d0
        stran = 0.d0
        if ((neps .eq. 6) .or. (neps .eq. 4)) then
            dstran(1:neps) = deps(1:neps)-BEHinteg%behavESVA%depsi_varc(1:neps)
            stran(1:neps) = epsm(1:neps)-BEHinteg%behavESVA%epsi_varc(1:neps)
        else if ((neps .eq. 3) .and. lCZM) then
! --------- No thermic strains for cohesive elements
            dstran(1:neps) = deps(1:neps)
            stran(1:neps) = epsm(1:neps)
        else if ((neps .eq. 12) .and. .not. lEpsa) then
! --------- For ENDO_HETEROGENE
            dstran(1:neps) = deps(1:neps)
            dstran(1:3) = dstran(1:3)-BEHinteg%behavESVA%depsi_varc(1:3)
        else if (lGradVari) then
! --------- For GRAD_VARI et GRAD_INCO
            ASSERT(neps .eq. 11 .or. neps .eq. 8)
            if (neps .eq. 11) then
                nepu = 6
            else if (neps .eq. 8) then
                nepu = 4
            end if
            dstran(1:nepu) = deps(1:nepu)-BEHinteg%behavESVA%depsi_varc(1:nepu)
            stran(1:nepu) = epsm(1:nepu)-BEHinteg%behavESVA%epsi_varc(1:nepu)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- epsm and deps become mechanical strains
        if (lGradVari) then
            epsm(1:nepu) = stran(1:nepu)
            deps(1:nepu) = dstran(1:nepu)
        else
            epsm(1:neps) = stran(1:neps)
            deps(1:neps) = dstran(1:neps)
        end if

        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare strains for integration: ', &
                neps, epsm(1:neps), deps(1:neps)
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! varcIsGEOM
!
! Detect 'GEOM' in external state variables
!
! Out lGeomInESVA      : flag for GEOM in external state variables (AFFE_VARC)
!
! --------------------------------------------------------------------------------------------------
    subroutine varcIsGEOM(lGeomInESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(inout) :: lGeomInESVA
! ----- Local
        character(len=8), parameter :: varc_geom = 'X'
        integer :: varc_indx
!   ------------------------------------------------------------------------------------------------
!
        lGeomInESVA = ASTER_FALSE

! ----- Detect 'GEOM' external state variables
        if (ca_nbcvrc_ .eq. 0) then
            lGeomInESVA = ASTER_FALSE
        else
            varc_indx = indik8(zk8(ca_jvcnom_), varc_geom, 1, ca_nbcvrc_)
            lGeomInESVA = varc_indx .ne. 0
        end if
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Detect GEOM in external state variables: ', lGeomInESVA
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isSolverIsExte
!
! Detect if external solver (MFront, UMAT) is used
!
! In  carcri           : parameters for comportment
! Out lMGIS            : logical for MFront/MGIS behaviour
! Out lUMAT            : logical for UMAT behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine isSolverIsExte(carcri, lMGIS, lUMAT)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        aster_logical, intent(out) :: lMGIS, lUMAT
!   ------------------------------------------------------------------------------------------------
        lMGIS = ASTER_FALSE
        lUMAT = ASTER_FALSE
        if (nint(carcri(EXTE_TYPE)) .eq. 1) then
!       MFront official
            lMGIS = ASTER_TRUE
        elseif (nint(carcri(EXTE_TYPE)) .eq. 2) then
!       MFront proto
            lMGIS = ASTER_TRUE
        elseif (nint(carcri(EXTE_TYPE)) .eq. 4) then
            lUMAT = ASTER_TRUE
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourOption
!
! Select objects to construct from option name
!
! In  option           : name of option to compute
! In  compor           : name of comportment definition (field)
! Out lMatr            : flag when tangent matrix
! Out lVect            : flag when internal forces vector
! Out lVari            : flag when internal state variables
! Out lSigm            : flag when stress and return code error
! Out codret           : return code when integrate behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourOption(option, compor, &
                               lMatr, lVect, &
                               lVari, lSigm, &
                               codret_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
        aster_logical, intent(out) :: lMatr, lVect, lVari, lSigm
        integer, optional, intent(out) :: codret_
! ----- Local
!   ------------------------------------------------------------------------------------------------
        aster_logical :: lPred
        integer :: copred, jv_copred, codret
!   ------------------------------------------------------------------------------------------------
        lVari = L_VARI(option)
        lSigm = L_SIGM(option)
        lVect = L_VECT(option)
        lMatr = L_MATR(option)
        lPred = L_PRED(option)
        if (lPred) then
            copred = 0
            if (compor(DEFO_LDC) .eq. 'MECANIQUE') then
                copred = 1
            end if
            call jevech('PCOPRED', 'E', jv_copred)
            zi(jv_copred) = copred
        end if
        codret = 0
        if (present(codret_)) then
            codret_ = codret
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module Behaviour_module
