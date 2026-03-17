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
!
subroutine exnode(mesh, zoneKeyword, nbUnilZone, nbNodeUnil, nolinoJv)
!
    implicit none
!
#include "asterfort/getnode.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: mesh
    character(len=16), intent(in) :: zoneKeyword
    integer(kind=8), intent(in) :: nbUnilZone
    character(len=24), intent(in) :: nolinoJv
    integer(kind=8), intent(in) :: nbNodeUnil
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get nodes from definition of LIAISON_UNILATER
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  nbUnilZone       : number of zones with LIAISON_UNILATER
! In  nbNodeUnil       : total number of nodes in LIAISON_UNILATER
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: listNode = '&&NBNODE.NOEU.NOEU'
    integer(kind=8) :: iUnilZone
    integer(kind=8), pointer :: nolino(:) => null(), nodeNume(:) => null()
    integer(kind=8) :: jdecal, nbNode, iNode
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Create object
    call wkvect(nolinoJv, 'V V I', nbNodeUnil, vi=nolino)

    jdecal = 1
    do iUnilZone = 1, nbUnilZone
        call getnode(mesh, zoneKeyword, iUnilZone, ' ', listNode, nbNode)
        if (nbNode .ne. 0) then
            call jeveuo(listNode, 'L', vi=nodeNume)
            do iNode = 1, nbNode
                nolino(jdecal) = nodeNume(iNode)
                jdecal = jdecal+1
            end do
        end if
        call jedetr(listNode)
    end do
!
    call jedema()
!
end subroutine
