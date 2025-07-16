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
subroutine te0097(option, nomte)
!
    implicit none
!
    character(len=16) :: option, nomte
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/elref2.h"
#include "asterfort/niinit.h"
#include "asterfort/sigvmc.h"
#include "asterfort/lteatt.h"
!
! --------------------------------------------------------------------------------------------------
!
!     BUT: CALCUL DES CONTRAINTES AUX POINTS DE GAUSS
!          ELEMENTS INCOMPRESSIBLE EN PETITES DEFORMATIONS
!
!          OPTION : 'SIEF_ELGA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgMax = 27, nbNodeMax = 27
    integer(kind=8) :: ndim, npg, jvWeightDisp, jvShapeDisp, jvDShapeDisp, jvShapePres
    integer(kind=8) :: idim
    integer(kind=8) :: jvSigm, jvDisp, jvGeom, jvMate, nbsig
    integer(kind=8) :: kpg, isig, iNodeDisp, iNodePres, iOSGS
    integer(kind=8) :: nbNodeDisp, nbNodePres, nbNodeGonf
    real(kind=8) :: sigmDisp(npgMax*6), angl_naut(3), instan, nharm, sigmTrac(npgMax)
    real(kind=8) :: presGaus(npgMax), dispU(3*nbNodeMax), dispP(nbNodeMax)
    integer(kind=8) :: vu(3, nbNodeMax), vg(nbNodeMax), vp(nbNodeMax), vpi(3, nbNodeMax)
    integer(kind=8) :: nbElrefe, iRefePres, iRefeGonf
    character(len=8) :: listElrefe(10), typmod(2)
    aster_logical :: lGonf
!
! --------------------------------------------------------------------------------------------------
!
    instan = r8vide()
    nharm = 0

! - Some unknonws
    if (lteatt('INCO', 'C2 ')) then
        lGonf = ASTER_FALSE
        iOSGS = 0
    elseif (lteatt('INCO', 'C2O')) then
        lGonf = ASTER_FALSE
        iOSGS = 1
    elseif (lteatt('INCO', 'C3 ')) then
        lGonf = ASTER_TRUE
        iOSGS = 0
    else
        ASSERT(ASTER_FALSE)
    end if

! - List of ELREFE
    call elref2(nomte, 10, listElrefe, nbElrefe)
    ASSERT(nbElrefe .ge. 2)
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)
    if (lGonf) then
        iRefePres = 3
        iRefeGonf = 2
    else
        iRefePres = 2
        iRefeGonf = 0
    end if

! - Get paramers of finite element for displacements
    call elrefe_info(elrefe=listElrefe(1), fami='RIGI', &
                     ndim=ndim, nno=nbNodeDisp, npg=npg, &
                     jpoids=jvWeightDisp, jvf=jvShapeDisp, jdfde=jvDShapeDisp)
    ASSERT(npg .le. npgMax)
    ASSERT(nbNodeDisp .le. nbNodeMax)

! - Get paramers of finite element for pres
    call elrefe_info(elrefe=listElrefe(iRefePres), fami='RIGI', &
                     nno=nbNodePres, &
                     jvf=jvShapePres)
    ASSERT(nbNodePres .le. nbNodeMax)

! - Get paramers of finite element for gonf
    nbNodeGonf = 0
    if (iRefeGonf .ne. 0) then
        call elrefe_info(elrefe=listElrefe(iRefeGonf), fami='RIGI', &
                         nno=nbNodeGonf)
    end if
    ASSERT(nbNodeGonf .le. nbNodeMax)

! - Modelling
    typmod = ' '
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get index of dof
    call niinit(typmod, &
                ndim, nbNodeDisp, nbNodeGonf, nbNodePres, iOSGS, &
                vu, vg, vp, vpi)

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMate)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Construct local anisotropic basis
    call getElemOrientation(ndim, nbNodeDisp, jvGeom, angl_naut)

! - Get displacements for u, v, w
    dispU = 0.d0
    do iNodeDisp = 1, nbNodeDisp
        do idim = 1, ndim
            dispU(idim+ndim*(iNodeDisp-1)) = zr(jvDisp-1+vu(iDim, iNodeDisp))
        end do
    end do

! - Get displacements for p
    dispP = 0.d0
    do iNodePres = 1, nbNodePres
        dispP(iNodePres) = zr(jvDisp-1+vp(iNodePres))
    end do

! - Compute stresses for displacement unknowns
    sigmDisp = 0.d0
    call sigvmc('RIGI', nbNodeDisp, ndim, nbsig, npg, &
                jvWeightDisp, jvShapeDisp, jvDShapeDisp, zr(jvGeom), dispU, &
                instan, angl_naut, zi(jvMate), nharm, sigmDisp)

    do kpg = 1, npg
        presGaus(kpg) = 0.d0
        do iNodePres = 1, nbNodePres
            presGaus(kpg) = presGaus(kpg)+ &
                            zr(jvShapePres-1+nbNodePres*(kpg-1)+iNodePres)*dispP(iNodePres)
        end do
        sigmTrac(kpg) = 0.d0
        do isig = 1, 3
            sigmTrac(kpg) = sigmTrac(kpg)+sigmDisp(nbsig*(kpg-1)+isig)
        end do
    end do

! - Output field
    call jevech('PCONTRR', 'E', jvSigm)

    do kpg = 1, npg
        do isig = 1, nbsig+1
            if (isig .le. 3) then
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmDisp(nbsig*(kpg-1)+isig)- &
                                                      sigmTrac(kpg)/3.d0+presGaus(kpg)
            elseif (isig .le. nbsig) then
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmDisp(nbsig*(kpg-1)+isig)
            else
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmTrac(kpg)/3.d0-presGaus(kpg)
            end if
        end do
    end do
!
end subroutine
