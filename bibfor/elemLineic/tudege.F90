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
subroutine tudege(lElga, nbNode, nbFourier, nbDof)
!
    use pipeElem_type
    use pipeElem_module
    use beamElem_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/ppgan2.h"
#include "jeveux.h"
!
    aster_logical, intent(in) :: lElga
    integer(kind=8), intent(in) :: nbNode, nbFourier, nbDof
!
! --------------------------------------------------------------------------------------------------
!
! Compute generalized strains for pipe element
!
! In  lElga            : flag for ELGA case
! In  nbNode           : number of nodes in element
! In  nbFourier        : number of Fourier modes
! In  nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvf, jdfde, jgano
    integer(kind=8) :: jvDispGlob
    integer(kind=8) :: jvDege
    integer(kind=8) :: iNode, iDof, iDege, nbDofByCell
    integer(kind=8) :: kpg, npg
    real(kind=8) :: dispLoca(nbDof)
    real(kind=8) :: degg(PIPE_NB_DEGE*PIPE_MAX_NPG), vectPg(PIPE_MAX_NPG), vectNo(PIPE_MAX_NODE)
    real(kind=8) :: hk, dhk
    type(pipeElem_Prop) :: pipeElem
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', npg=npg, &
                     jvf=jvf, jdfde=jdfde, jgano=jgano)
    ASSERT(npg .le. PIPE_MAX_NPG)
    nbDofByCell = (6+3+6*(nbFourier-1))

! - Output field
    if (lElga) then
        call jevech('PDEFOPG', 'E', jvDege)
    else
        call jevech('PDEFOGR', 'E', jvDege)
    end if

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)
    if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
        pipeElem%beamElem%elemLength = pipeElem%beamElem%thetaElbow*pipeElem%beamElem%radiusElbow
    end if

! - Get displacements in local base
    call jevech('PDEPLAR', 'L', jvDispGlob)
    call pipeGetDisp(pipeElem, nbDof, zr(jvDispGlob), dispLoca)

! - Loop on Gauss points (on segment)
    degg = 0.d0
    do kpg = 1, npg
        do iNode = 1, nbNode
! --------- Shape function and its derivative
            hk = zr(jvf-1+nbNode*(kpg-1)+iNode)
            dhk = zr(jdfde-1+nbNode*(kpg-1)+iNode)*(2.d0/pipeElem%beamElem%elemLength)

! --------- Compute DEGE_ELGA
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_EPXX) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_EPXX)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DX)
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAXY) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAXY)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DY)- &
                hk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DRZ)
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAXZ) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAXZ)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DZ)+ &
                hk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DRY)
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAT) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_GAT)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DRX)
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_KY) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_KY)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DRY)
            degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_KZ) = &
                degg(PIPE_NB_DEGE*(kpg-1)+PIPE_DEGE_KZ)+ &
                dhk*dispLoca(nbDofByCell*(iNode-1)+PIPE_DOF_DRZ)
        end do
    end do

    do iDege = 1, PIPE_NB_DEGE
        do kpg = 1, npg
            vectPg(kpg) = degg(PIPE_NB_DEGE*(kpg-1)+iDege)
            if (lElga) then
                zr(jvDege+PIPE_NB_DEGE*(kpg-1)+iDege-1) = vectPg(kpg)
            end if
        end do
        if (.not. lElga) then
            call ppgan2(jgano, 1, 1, vectPg, vectNo)
            do iNode = 1, nbNode
                zr(jvDege+PIPE_NB_DEGE*(iNode-1)+iDege-1) = vectNo(iNode)
            end do
        end if
    end do
!
end subroutine
