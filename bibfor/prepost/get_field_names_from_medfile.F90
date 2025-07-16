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
subroutine get_field_names_from_medfile(idfimd, field_name_vector)
! person_in_charge: nicolas.sellenet at edf.fr
!-----------------------------------------------------------------------
!    Read all field names in a med file
!
!     In:
!        idfimd : med file id
!
!     In/Out:
!        field_name_vector : name of jeveux vector to fill
!_____________________________________________________________________
!
    use as_med_module, only: as_med_open
    implicit none
!
! 0.1. ==> Arguments
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_mfdfdi.h"
#include "asterfort/as_mfdnfc.h"
#include "asterfort/as_mfdnfd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    med_idt idfimd
    character(len=*) :: field_name_vector
!
    integer(kind=8) :: nbcham, nbcmp, jcmp
    integer(kind=8) :: tycha, junit
    integer(kind=8) :: nseqca
    integer(kind=8) :: iret, i
    character(len=80), pointer :: v_names(:) => null()
!
    character(len=64) :: nomcha
!
    call jemarq()
!
    call jedetr(field_name_vector)
!
    call as_mfdnfd(idfimd, nbcham, iret)
    call wkvect(field_name_vector, 'V V K80', nbcham, vk80=v_names)
    do i = 1, nbcham
        call as_mfdnfc(idfimd, i, nbcmp, iret)
        call wkvect('&&LRCEME.NOMCMP_K16', 'V V K16', nbcmp, jcmp)
        call wkvect('&&LRCEME.UNITCMP', 'V V K16', nbcmp, junit)
        call as_mfdfdi(idfimd, i, nomcha, tycha, zk16(jcmp), &
                       zk16(junit), nseqca, iret)
        v_names(i) = nomcha
        call jedetr('&&LRCEME.NOMCMP_K16')
        call jedetr('&&LRCEME.UNITCMP')
    end do
!
    call jedema()
!
end subroutine
