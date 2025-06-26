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

subroutine cnsreddepl(cns1z)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cnsred.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*), intent(inout) :: cns1z
! ---------------------------------------------------------------------
! BUT: REDUIRE le champ à ces composantes concernées par le calcul du
!      travail exterieur : DX, DY, DZ, DRX, DRY et DRZ
! ---------------------------------------------------------------------
!     ARGUMENTS:
! CNS1Z  IN/OUT  K19 : SD CHAM_NO_S
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: k, ncmp, ncmp1, j
    character(len=19) :: cns1
    aster_logical :: l_cmp_present(6)
    character(len=8), pointer :: cnsc1(:) => null()
    integer(kind=8), pointer :: cnsd1(:) => null()
    character(len=8) :: depl_cmp(6), compSelect(6)
    data depl_cmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
!     ------------------------------------------------------------------
    call jemarq()
!
    cns1 = cns1z
!
    call jeveuo(cns1//'.CNSD', 'L', vi=cnsd1)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cnsc1)
    ncmp1 = cnsd1(2)
!
    l_cmp_present(:) = .false.
    do k = 1, ncmp1
        do j = 1, 6
            if (cnsc1(k) .eq. depl_cmp(j)) then
                l_cmp_present(j) = .true.
            end if
        end do
    end do
!
    ncmp = 0
    do j = 1, 6
        if (l_cmp_present(j)) then
            ncmp = ncmp+1
            compSelect(ncmp) = depl_cmp(j)
        end if
    end do
!
    call cnsred(cns1, 0, [0], ncmp, compSelect, &
                'V', cns1)

    call jedema()
end subroutine
