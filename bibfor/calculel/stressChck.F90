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
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.a
! --------------------------------------------------------------------
!
subroutine stressChck(stressRefeZ, stressZ, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cestas.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: stressRefeZ, stressZ
    integer, intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Check the sous-points of the input stress field
!
! Is stress field is correct for the number of sous points ?
!
! --------------------------------------------------------------------------------------------------
!
! In  stressRefe       : reference for stress field
! In  stress           : stress field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: stressS = '&&SGCOMP.SIGM_S'
    character(len=19), parameter :: stressRefeS = '&&SGCOMP.SIGR_S'
    character(len=8) :: meshRefe, mesh
    aster_logical :: noSameNpg, noSameCmp
    integer :: iad1, iad2
    integer :: nbCell, iCell, vali(3)
    integer :: npgRefe, npg
    integer :: ncmpRefe, ncmp
    character(len=8) :: cellName
    integer :: jvCesdRefe, jvCeslRefe
    integer :: jvCesd, jvCesl
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    iret = 0

! - Check meshes
    call dismoi('NOM_MAILLA', stressRefeZ, 'CHAMP', repk=meshRefe)
    call dismoi('NOM_MAILLA', stressZ, 'CHAMP', repk=mesh)
    if (meshRefe .ne. mesh) then
        call utmess('F', 'FIELD0_10')
    end if

! - Create reduced fields for stress
    call celces(stressZ, 'V', stressS)
    call cestas(stressS)
    call celces(stressRefeZ, 'V', stressRefeS)
    call cestas(stressRefeS)

! - Access to fields
    call jeveuo(stressRefeS//'.CESD', 'L', jvCesdRefe)
    call jeveuo(stressRefeS//'.CESL', 'L', jvCeslRefe)
    call jeveuo(stressS//'.CESD', 'L', jvCesd)
    call jeveuo(stressS//'.CESL', 'L', jvCesl)

    nbCell = zi(jvCesd-1+2)

! - Check
    noSameNpg = ASTER_FALSE
    noSameCmp = ASTER_FALSE
    do iCell = 1, nbCell
        call cesexi('C', jvCesdRefe, jvCeslRefe, iCell, 1, 1, 1, iad1)
        if (iad1 .le. 0) then
            cycle
        end if
        npgRefe = zi(jvCesdRefe-1+5+4*(iCell-1)+1)
        ncmpRefe = zi(jvCesdRefe-1+5+4*(iCell-1)+3)

        call cesexi('C', jvCesd, jvCesl, iCell, 1, 1, 1, iad2)
        if (iad2 .gt. 0) then
            npg = zi(jvCesd-1+5+4*(iCell-1)+1)
            ncmp = zi(jvCesd-1+5+4*(iCell-1)+3)
! --------- Check number of Gauss points
            if (npgRefe .ne. 0 .and. npg .ne. 0) then
                if (npgRefe .ne. npg) then
                    call jenuno(jexnum(mesh//'.NOMMAI', iCell), cellName)
                    vali(1) = npgRefe
                    vali(2) = npg
                    call utmess('I', 'FIELD0_11', sk=cellName, ni=2, vali=vali)
                    noSameNpg = ASTER_TRUE
                    goto 40
                end if
            end if
! --------- Check number of components
            if (ncmpRefe .ne. 0 .and. ncmp .ne. 0) then
                if (ncmpRefe .ne. ncmp) then
                    call jenuno(jexnum(mesh//'.NOMMAI', iCell), cellName)
                    vali(1) = ncmpRefe
                    vali(2) = ncmp
                    call utmess('I', 'FIELD0_12', sk=cellName, ni=2, vali=vali)
                    noSameCmp = ASTER_TRUE
                    goto 40
                end if
            end if
        end if
40      continue
    end do
!
    if (noSameNpg .or. noSameCmp) then
        iret = 1
        call utmess("E", 'FIELD0_13')
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', stressS)
    call detrsd('CHAM_ELEM_S', stressRefeS)
!
    call jedema()
!
end subroutine
