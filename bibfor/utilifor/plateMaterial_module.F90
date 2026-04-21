! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! aslint: disable=W0413
!
! ==================================================================================================
!
! Module for management of material of plates
!
! ==================================================================================================
!
module plateMaterial_module
! ==================================================================================================
    use MaterialPara_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: getMultiLayerNbLayer
    private :: chckMultiLayer
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! chckMultiLayer
!
! Check properties of multi-layered plates
!
! --------------------------------------------------------------------------------------------------
    subroutine chckMultiLayer(materPara, nbLayer)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
        integer(kind=8), intent(in) :: nbLayer
! ----- Local
        character(len=16), parameter :: elasKeyword = 'ELAS_COQMU'
        character(len=8), parameter :: fami = 'FPG1'
        character(len=1), parameter :: poum = '+'
        integer(kind=8), parameter :: kpg = 1, ksp = 1
        real(kind=8), parameter :: r8bid = 0.d0
        integer(kind=8) :: iLayer, jvCacoqu, iret, jvNbspIn
        real(kind=8) :: epTotaSum, epTota, epLayer(1)
        character(len=3) :: iLayerStr
        character(len=2) :: iValeStr
        character(len=16) :: propName
        integer(kind=8) :: propCode(1)
!   ------------------------------------------------------------------------------------------------
!
        call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jvNbspIn)
        ASSERT(iret .eq. 0)
        iLayer = 0
        epTotaSum = 0.d0
        epLayer(1) = 0.d0
        call jevech('PCACOQU', 'L', jvCacoqu)
        epTota = zr(jvCacoqu)
5       continue
        iLayer = iLayer+1
        call codent(iLayer, 'G', iLayerStr)
        call codent(1, 'G', iValeStr)
        propName = 'C'//iLayerStr//'_V'//iValeStr
        call rcvalb(fami, kpg, ksp, poum, &
                    materPara%jvMaterCode, ' ', elasKeyword, &
                    0, ' ', [r8bid], &
                    1, propName, epLayer, propCode(1), 0)
        if (propCode(1) .eq. 0) then
            epTotaSum = epTotaSum+epLayer(1)
            goto 5
        end if
        if (epTotaSum .ne. 0.d0) then
            if ((iLayer-1) .ne. nbLayer) then
                call utmess('F', 'PLATE1_51', ni=2, vali=[iLayer-1, nbLayer])
            end if
            if (abs(epTota-epTotaSum)/epTota .gt. 1.d-2) then
                call utmess('F', 'PLATE1_52', nr=2, valr=[epTotaSum, epTota])
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMultiLayerNbLayer
!
! Get number of layers
!
! --------------------------------------------------------------------------------------------------
    subroutine getMultiLayerNbLayer(materPara, lDKTG, nbLayer)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
        aster_logical, intent(in) :: lDKTG
        integer(kind=8), intent(out) :: nbLayer
! ----- Local
        integer(kind=8) :: jvNbsp
!   ------------------------------------------------------------------------------------------------
!
        nbLayer = 0
        if (lDKTG) then
            nbLayer = 1
        else
            call jevech('PNBSP_I', 'L', jvNbsp)
            nbLayer = zi(jvNbsp-1+1)
            if (nbLayer .le. 0) then
                call utmess('F', 'ELEMENTS_46')
            end if
        end if
        call chckMultiLayer(materPara, nbLayer)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module plateMaterial_module
