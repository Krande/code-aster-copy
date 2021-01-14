! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
module SolidShell_NonLinear_Hexa_module
! ==================================================================================================
use SolidShell_type
use SolidShell_Debug_module
use SolidShell_Geometry_Hexa_module
use SolidShell_Kinematic_Hexa_module
use SolidShell_Stabilization_Hexa_module
use SolidShell_Utilities_module
use Behaviour_module
! ==================================================================================================
implicit none
! ==================================================================================================
public  :: compNonLinearHexa
private :: compSmallStrainHexa
! ==================================================================================================
private
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/SolidShell_type.h"
#include "asterc/isnnem.h"
#include "asterfort/assert.h"
#include "asterfort/btsig.h"
#include "asterfort/codere.h"
#include "asterfort/jevech.h"
#include "asterfort/nmcomp.h"
#include "asterfort/tecach.h"
#include "blas/dcopy.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! compNonLinearHexa
!
! Compute non-linear options for HEXA
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  matePara         : parameters of material
! In  behaPara         : parameters of behaviour
!
! --------------------------------------------------------------------------------------------------
subroutine compNonLinearHexa(option, elemProp, cellGeom, matePara, behaPara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    character(len=16), intent(in)   :: option
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    type(SSH_MATE_PARA), intent(in) :: matePara
    type(SSH_BEHA_PARA), intent(in) :: behaPara
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    integer :: nbIntePoint, nbVari
    integer :: jtab(7), iret
    integer :: jvGeom, jvMater
    integer :: jvMatr, jvVect, jvSigmP, jvVariP, jvVariX, jvCodret
    integer :: jvTimeM, jvTimeP, jvSigmM, jvVariM, jvDispM, jvDispIncr
!   ------------------------------------------------------------------------------------------------
!

! - Properties of finite element
    nbIntePoint = elemProp%elemInte%nbIntePoint

! - Get input fields
    jvGeom  = cellGeom%jvGeom
    jvMater = matePara%jvMater
    call jevech('PINSTMR', 'L', jvTimeM)
    call jevech('PINSTPR', 'L', jvTimeP)
    call jevech('PCONTMR', 'L', jvSigmM)
    call jevech('PDEPLMR', 'L', jvDispM)
    call jevech('PDEPLPR', 'L', jvDispIncr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    jvVariM = jtab(1)
    nbVari  = max(jtab(6), 1)*jtab(7)

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Get output fields
    jvMatr = isnnem()
    if (behaPara%lMatr) then
        if (behaPara%lMatrSyme) then
            call jevech('PMATUUR', 'E', jvMatr)
        else
! --------- Matrices non symétriques pas prises en charge
            ASSERT(ASTER_FALSE)
            !call jevech('PMATUNS', 'E', jvMatr)
        endif
    endif
    jvVect = isnnem()
    if (behaPara%lVect) then
        call jevech('PVECTUR', 'E', jvVect)
    endif
    jvSigmP = isnnem()
    if (behaPara%lSigm) then
        call jevech('PCONTPR', 'E', jvSigmP)
        call jevech('PCODRET', 'E', jvCodret)
    endif
    jvVariP = isnnem()
    if (behaPara%lVari) then
        call jevech('PVARIPR', 'E', jvVariP)
        call jevech('PVARIMP', 'L', jvVariX)
        call dcopy(nbIntePoint*nbVari, zr(jvVariX), 1, zr(jvVariP), 1)
    endif

! - Compute
    if (behaPara%defoComp .eq. 'PETIT') then
        call compSmallStrainHexa(option      , elemProp      , cellGeom, geomHexa,&
                                 matePara    , behaPara      ,&
                                 nbIntePoint , nbVari        ,&
                                 zr(jvTimeM) , zr(jvTimeP)   ,&
                                 zr(jvDispM) , zr(jvDispIncr),&
                                 zr(jvSigmM) , zr(jvVariM)   ,&
                                 zr(jvSigmP) , zr(jvVariP)   ,&
                                 zr(jvMatr)  , zr(jvVect)    ,&
                                 zi(jvCodret))

    else
        ASSERT(ASTER_FALSE)
    endif
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSmallStrainHexa
!
! Compute non-linear options for HEXA - Small strains
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  geomHexa         : geometric properties for HEXA cell
! In  matePara         : parameters of material
! In  behaPara         : parameters of behaviour
! In  nbIntePoint      : number of integration points on cell
! In  nbVari           : number of internal state variables
! In  timePrev         : time of previous time step
! In  timeCurr         : time of current time step
! In  dispPrev         : displacement at beginning of current time step
! In  dispIncr         : displacement since beginning of current time step
! In  sigm             : stress at beginning of current time step (before integration)
! In  vim              : internal state vari. at beginning of current time step (before integration)
! Out sigp             : stress at end of current time step (after integration)
! Out vip              : internal state vari. at end of current time step (after integration)
! Out matr             : tangent matrix
! Out vect             : vector of internal forces
! Out codret           : error code from integration of behaviour
!
! --------------------------------------------------------------------------------------------------
subroutine compSmallStrainHexa(option     , elemProp, cellGeom, geomHexa,&
                               matePara   , behaPara,&
                               nbIntePoint, nbVari  ,&
                               timePrev   , timeCurr,&
                               dispPrev   , dispIncr,&
                               sigm       , vim     ,&
                               sigp       , vip     ,&
                               matr       , vect    ,&
                               codret)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    character(len=16), intent(in)   :: option
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    type(SSH_GEOM_HEXA), intent(in) :: geomHexa
    type(SSH_MATE_PARA), intent(in) :: matePara
    type(SSH_BEHA_PARA), intent(in) :: behaPara
    integer, intent(in)      :: nbIntePoint, nbVari
    real(kind=8), intent(in) :: timePrev, timeCurr
    real(kind=8), intent(in) :: dispPrev(SSH_NBDOF_HEXA), dispIncr(SSH_NBDOF_HEXA)
    real(kind=8), intent(in) :: sigm(SSH_SIZE_TENS, nbIntePoint), vim(nbVari, nbIntePoint)
    real(kind=8), intent(out) :: sigp(SSH_SIZE_TENS, nbIntePoint), vip(nbVari, nbIntePoint)
    real(kind=8), intent(out) :: matr(SSH_SIZE_MATR_HEXA), vect(SSH_NBDOF_HEXA)
    integer, intent(out) :: codret
! - Local
    integer, parameter :: ksp = 1
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer :: cod(nbIntePoint), kpg, iTens, nbDof, nbDofGeom
    integer :: jvCoor, jvWeight
    type(SSH_KINE_HEXA) :: kineHexa
    type(SSH_STAB_HEXA) :: stabHexa
    real(kind=8) :: zeta, poids, jacob
    real(kind=8) :: UeffKpg, Ueff
    real(kind=8) :: geomCurr(SSH_NBDOFG_HEXA), dispCurr(SSH_NBDOF_HEXA)
    real(kind=8) :: epsiPrev(SSH_SIZE_TENS), epsiIncr(SSH_SIZE_TENS)
    real(kind=8) :: sigmPost(SSH_SIZE_TENS), sigmPrep(SSH_SIZE_TENS)
    real(kind=8), dimension(SSH_NBDOF_HEXA, SSH_NBDOF_HEXA) :: matrMate
    real(kind=8), dimension(SSH_SIZE_TENS, SSH_SIZE_TENS) :: dsidep
!   ------------------------------------------------------------------------------------------------
!
    if (behaPara%lSigm) then
        codret = 0
        sigp   = 0.d0
    endif
    if (behaPara%lVect) then
        vect = 0.d0
    endif
    if (behaPara%lMatr) then
        ASSERT(behaPara%lMatrSyme)
        matr = 0.d0
    endif

! - Properties of finite element
    nbDof     = elemProp%nbDof
    nbDofGeom = elemProp%nbDofGeom
    jvCoor    = elemProp%elemInte%jvCoor
    jvWeight  = elemProp%elemInte%jvWeight

! - Properties of finite element
    ASSERT(.not.behaPara%lLarge)

! - Total displacement from initial configuration
    dispCurr = dispIncr + dispPrev

! - Update configuration: initial
    geomCurr(1:SSH_NBDOFG_HEXA) = cellGeom%geomInit(1:SSH_NBDOFG_HEXA)

! - Compute gradient matrix in covariant basis
    kineHexa%lLarge = ASTER_FALSE
    call compBCovaMatrHexa(geomCurr, kineHexa)

! - Compute gradient matrix in cartesian frame
    call compBCartMatrHexa(geomHexa, kineHexa)
    if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_ = ASTER_TRUE)

! - Loop on Gauss points
    Ueff = 0.d0
    do kpg = 1, nbIntePoint
        zeta  = zr(jvCoor-1+3*kpg)
        poids = zr(jvWeight-1+kpg)
        jacob = poids * cellGeom%detJac0

! ----- Compute EAS B matrix in cartesian frame at current Gauss point
        call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! ----- Compute B matrix
        call compBMatrHexa(zeta, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_ = ASTER_TRUE)

! ----- Compute small strains at beginning of time step
        call compEpsiHexa(kineHexa, dispPrev, epsiPrev)

! ----- Compute increment of small strains
        call compEpsiHexa(kineHexa, dispIncr, epsiIncr)

! ----- Pre-treatment of stresses and strains
        epsiPrev(4) = epsiPrev(4)/rac2
        epsiPrev(5) = epsiPrev(5)/rac2
        epsiPrev(6) = epsiPrev(6)/rac2
        epsiIncr(4) = epsiIncr(4)/rac2
        epsiIncr(5) = epsiIncr(5)/rac2
        epsiIncr(6) = epsiIncr(6)/rac2
        do iTens = 1, 3
            sigmPrep(iTens)   = sigm(iTens, kpg)
            sigmPrep(iTens+3) = sigm(iTens+3, kpg)*rac2
        end do

! ----- Integrate behaviour law
        sigmPost = 0.d0
        dsidep   = 0.d0
        cod(kpg) = 0
        call nmcomp(behaPara%BEHInteg,&
                    elemProp%elemInte%inteFami, kpg, ksp,&
                    SSH_NDIM, typmod,&
                    matePara%jvMater,&
                    zk16(behaPara%jvCompor),&
                    zr(behaPara%jvCarcri),&
                    timePrev, timeCurr,&
                    SSH_SIZE_TENS, epsiPrev, epsiIncr,&
                    SSH_SIZE_TENS, sigmPrep,&
                    vim(1, kpg), option, matePara%mateBase,&
                    sigmPost, vip(1, kpg),&
                    SSH_SIZE_TENS*SSH_SIZE_TENS, dsidep,&
                    cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 99
        endif

! ----- Post-treatment of stresses and matrix
        if (behaPara%lSigm) then
            sigmPost(4)     = sigmPost(4)/rac2
            sigmPost(5)     = sigmPost(5)/rac2
            sigmPost(6)     = sigmPost(6)/rac2
        endif
        if (behaPara%lMatr) then
            dsidep(4:6, 4:6) = dsidep(4:6, 4:6)/2.d0
            dsidep(4:6, 1:3) = dsidep(4:6, 1:3)/rac2
            dsidep(1:3, 4:6) = dsidep(1:3, 4:6)/rac2
        endif

! ----- Compute effective shear modulus for stabilization
        call compStabModulusHexa(sigmPost, epsiIncr, dsidep, UeffKpg)
        Ueff = Ueff + UeffKpg*poids/8.d0

! ----- Compute material part  at current Gauss point
        if (behaPara%lMatr) then
            matrMate = 0.d0
            call prodBTDB(dsidep, SSH_SIZE_TENS, elemProp%nbDof, kineHexa%B, matrMate)
        endif

! ----- Update tangent matrix 
        if (behaPara%lMatr) then
            ASSERT(behaPara%lMatrSyme)
            call updateMatrSyme(nbDof, matrMate, matr, jacob)
        endif

! ----- Update internal force at current Gauss point
        if (behaPara%lVect) then
            call btsig(elemProp%nbDof, SSH_SIZE_TENS, jacob,&
                       kineHexa%B    , sigmPost     , vect)
        endif

! ----- Save stresses
        if (behaPara%lSigm) then
            do iTens = 1, SSH_SIZE_TENS
                sigp(iTens, kpg) = sigmPost(iTens)
            end do
        endif

    end do

! - Stabilization
    if (behaPara%lVect) then
        call compStabSigmHexa(geomHexa, kineHexa, Ueff, dispCurr, stabHexa)
    endif
    if (behaPara%lMatr) then
        call compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
    endif
    if (behaPara%lVect) then
        call compStabForcHexa(kineHexa, stabHexa)
    endif

! - Save matrix and vector
    if (behaPara%lMatr) then
        call updateMatrSyme(nbDofGeom, stabHexa%matrStabMate, matr)
    endif
    if (behaPara%lVect) then
        vect(1:SSH_NBDOFG_HEXA) =&
            vect(1:SSH_NBDOFG_HEXA) + cellGeom%detJac0*stabHexa%forcStab
    endif
!
 99 continue

! - Return code summary
    if (behaPara%lSigm) then
        call codere(cod, nbIntePoint, codret)
    endif
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_NonLinear_Hexa_module
