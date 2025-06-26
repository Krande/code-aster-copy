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
subroutine asseVectField(vectAsse, numeDof, vectScalType, &
                         nbVectElem, listVectElem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: vectAsse
    character(len=14), intent(in) :: numeDof
    integer(kind=8), intent(in) :: vectScalType
    integer(kind=8), intent(in) :: nbVectElem
    character(len=*), intent(in) :: listVectElem(nbVectElem)
!
! --------------------------------------------------------------------------------------------------
!
! Assembly vector when elementary vector is a nodal field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: vectAsseWork = '&&ASSVEC.CHAMNO'
    integer(kind=8) :: iVectElem, iNode, iEqua
    integer(kind=8) :: iexi, iret
    integer(kind=8) :: nbNode, nbEqua
    character(len=19) :: vectElem, nodeField
    character(len=24), pointer :: relr(:) => null()
    character(len=1) :: ktyp
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: valeWork(:) => null()

!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    do iVectElem = 1, nbVectElem
        vectElem = listVectElem(iVectElem)
        call jeexin(vectElem//'.RELR', iexi)
        if (iexi .eq. 0) cycle
        call jeveuo(vectElem//'.RELR', 'L', vk24=relr)
        call jelira(vectElem//'.RELR', 'LONUTI', nbNode)
        do iNode = 1, nbNode
            nodeField = relr(iNode) (1:19)
            call jeexin(nodeField//'.VALE', iexi)
            if (iexi .gt. 0) then
                call jeveuo(vectAsse//'.VALE', 'E', vr=vale)
                call jelira(vectAsse//'.VALE', 'TYPE', cval=ktyp)
                ASSERT(ktyp .eq. 'R')
                ASSERT(vectScalType .eq. 1)
                call vtcreb(vectAsseWork, 'V', ktyp, &
                            nume_ddlz=numeDof, &
                            nb_equa_outz=nbEqua)
                call vtcopy(nodeField, vectAsseWork, iret)
                if (iret .ne. 0) then
                    call utmess("F", "FIELD0_5")
                end if
                call jeveuo(vectAsseWork//'.VALE', 'L', vr=valeWork)
                do iEqua = 1, nbEqua
                    vale(iEqua) = vale(iEqua)+valeWork(iEqua)
                end do
                call detrsd('CHAM_NO', vectAsseWork)
            end if
        end do
    end do
!
    call jedema()
end subroutine
