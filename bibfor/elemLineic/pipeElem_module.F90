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
! ==================================================================================================
!
! Module for pipe elements
!
! ==================================================================================================
!
module pipeElem_module
! ==================================================================================================
    use beamElem_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private :: pipeGetType, pipeGetSection
    public :: pipeCheckMetric
! ==================================================================================================
    private
#include "asterc/r8pi.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/poutre_modloc.h"
#include "jeveux.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! pipeGetSection
!
! Get section of pipe
!
! Out sectPipe         : properties of pipe's section
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetSection(sectPipe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(sectPipe_Prop), intent(out) :: sectPipe
! ----- Local
        integer(kind=8), parameter :: nbCara = 2
        character(len=8), parameter :: caraName(nbCara) = (/'R1 ', 'EP1'/)
        real(kind=8) :: caraVale(nbCara)
!   ------------------------------------------------------------------------------------------------
!
        call poutre_modloc('CAGEP1', caraName, nbCara, lvaleur=caraVale)
        sectPipe%radiusExt = caraVale(1)
        sectPipe%thickness = caraVale(2)
        sectPipe%radiusMoy = sectPipe%radiusExt-sectPipe%thickness/2.d0
        sectPipe%radiusInt = sectPipe%radiusExt-sectPipe%thickness
        sectPipe%area = r8pi()*(sectPipe%radiusExt**2-sectPipe%radiusInt**2)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetType
!
! Get type of pipe
!
! In  jvPipe           : adress for parameters of pipe
! In  nbNode           : number of nodes of element
! Out pipeType         : type of pipe
! Out lModiMetric      : using MODI_METRIQUE
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetType(jvPipe, nbNode, pipeType, lModiMetric)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvPipe
        integer(kind=8), intent(in) :: nbNode
        integer(kind=8), intent(out) :: pipeType
        aster_logical, intent(out) :: lModiMetric
! ----- Local
        integer(kind=8) :: iCoude, iCmp
!   ------------------------------------------------------------------------------------------------
!
        pipeType = PIPE_TYPE_UNDEF
        if (nbNode .eq. 3) then
            iCmp = 9
        else if (nbNode .eq. 4) then
            iCmp = 12
        else
            ASSERT(ASTER_FALSE)
        end if
        iCoude = nint(zr(jvPipe-1+iCmp+1))
        lModiMetric = iCoude .lt. 10
        if (iCoude .ge. 10) then
            iCoude = iCoude-10
        end if
        if (iCoude .eq. 0) then
            pipeType = PIPE_TYPE_STRAIGHT
        elseif (iCoude .eq. 1) then
            pipeType = PIPE_TYPE_ELBOW
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeCheckMetric
!
! Check properties of pipe for modiMetric
!
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeCheckMetric(nomteZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: nomteZ
! ----- Local
        integer(kind=8) :: jvPipe, jvCodret, jvCheck, jvIndicR
        integer(kind=8) :: nbFourier, nbNode
        integer(kind=8) :: pipeType
        real(kind=8) :: criteria, tole
        type(sectPipe_Prop) :: sectPipe
        aster_logical :: lModiMetric
!   ------------------------------------------------------------------------------------------------
!
        call elrefe_info(fami='RIGI', nno=nbNode)
        call jevech('PCAORIE', 'L', jvPipe)
        call jevech('PCHCKPR', 'L', jvCheck)
        call jevech('PCODRET', "E", jvCodret)
        call jevech('PINDICR', "E", jvIndicR)
        nbFourier = 3
        if (nomteZ .eq. 'MET6SEG3') then
            nbFourier = 6
        end if
        call pipeGetType(jvPipe, nbNode, pipeType, lModiMetric)
        call pipeGetSection(sectPipe)
        tole = zr(jvCheck-1+1)
        criteria = sectPipe%thickness/sectPipe%radiusMoy
        zr(jvIndicR-1+1) = criteria
        if (lModiMetric) then
            if (criteria .gt. tole) then
                zi(jvCodret-1+1) = 2
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module pipeElem_module
