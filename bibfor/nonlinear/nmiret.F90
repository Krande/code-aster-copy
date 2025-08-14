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
subroutine nmiret(codretZ, tabret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmpi_comm_logical.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdmpic.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: codretZ
    aster_logical, intent(out) :: tabret(0:10)
!
! --------------------------------------------------------------------------------------------------
!
! Treatement of return code for elementary test
!
! --------------------------------------------------------------------------------------------------
!
! In  codretZ     : field for return code
! Out tabret      : flags
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: chamns = '&&NMIRET.CHAMNS'
    integer(kind=8) :: iret, jcesd, jcesl, nbCell, nbCmp
    integer(kind=8) :: iCell, iad
    character(len=8) :: physQuan, mesh
    character(len=19) :: codret
    integer(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: cesk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    codret = codretZ
    tabret = ASTER_FALSE

! - Check input field
    call exisd('CHAM_ELEM', codret, iret)
    if (iret .ne. 0) then
        call sdmpic('CHAM_ELEM', codret)

! ----- Conversion and access
        call celces(codret, 'V', chamns)
        call jeveuo(chamns(1:19)//'.CESK', 'L', vk8=cesk)
        call jeveuo(chamns(1:19)//'.CESD', 'L', jcesd)
        call jeveuo(chamns(1:19)//'.CESV', 'L', vi=cesv)
        call jeveuo(chamns(1:19)//'.CESL', 'L', jcesl)

! ----- Get parameters
        ASSERT(zi(jcesd-1+3) .eq. 1)
        ASSERT(zi(jcesd-1+4) .eq. 1)
        physQuan = cesk(2)
        ASSERT(physQuan .eq. 'CODE_I')
        nbCell = zi(jcesd-1+1)
        nbCmp = zi(jcesd-1+2)
        ASSERT(nbCmp .eq. 1)

! ----- Set flags
        do iCell = 1, nbCell
            call cesexi('C', jcesd, jcesl, iCell, 1, 1, nbCmp, iad)
            if (iad .le. 0) cycle
            iret = cesv(iad)
            if (iret .eq. 0) then
            else if (iret .lt. 11 .and. iret .gt. 0) then
                tabret(iret) = ASTER_TRUE
            else
                ASSERT(ASTER_FALSE)
            end if
        end do

! ----- Il faut faire une synthèse en HPC
        call dismoi('NOM_MAILLA', codret, 'CHAM_ELEM', repk=mesh)
        if (isParallelMesh(mesh)) then
            call asmpi_comm_logical("MPI_LOR", nbval=10, vl=tabret)
        end if
        do iret = 1, 10
            if (tabret(iret)) then
                tabret(0) = ASTER_TRUE
            end if
        end do
!
        call detrsd('CHAM_ELEM_S', chamns)
    end if
!
    call jedema()
end subroutine
