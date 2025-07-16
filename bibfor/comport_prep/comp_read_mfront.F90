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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_read_mfront(keywf, i_comp, extern_addr)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/jeveuo.h"
!
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: i_comp
    character(len=16), intent(out) :: extern_addr
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get parameters for external programs (MFRONT)
!
! --------------------------------------------------------------------------------------------------
!
! In      keywf      : factor keyword to read (COMPORTEMENT)
! In      i_comp     : factor keyword index
! Out     extern_ptr : MGIS address
!
! --------------------------------------------------------------------------------------------------
    character(len=8) :: mgb
    character(len=16), pointer :: addr(:) => null()
    integer(kind=8) :: nbret
!
! - Get parameters
!
    ASSERT(i_comp .ne. 0)

    call getvid(keywf, "COMPOR_MFRONT", i_comp, scal=mgb, nbret=nbret)
    if (nbret .ne. 1) then
        ! May happen when called from comp_ntvari if iMap does not match a factkeyword!
        print *, "MGISDBG: compor_mgis object not found: factor keyword: '", &
            keywf, "', occ:", i_comp
        ! hmm...
        ASSERT(.false.)
    end if

    call jeveuo(mgb//'.ADDR', 'L', vk16=addr)
    extern_addr = addr(1)
!
end subroutine
