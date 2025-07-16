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
subroutine afva01(affeType, dsName, fieldType, lautr)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cmpcha.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: affeType, dsName
    character(len=16), intent(in) :: fieldType
    aster_logical, intent(out) :: lautr
!
! ----------------------------------------------------------------------
! BUT : DIRE SI DANS LA SD NOMSD DE TYPE TYPSD=CHAMP/EVOL+NOMSYM
!       ON TROUVE DES COMPOSANTES AUTRES QUE 'TEMP' ET 'LAGR'
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nb_cmp, k, jordr, j, iret, nbordr(1), ibid
    character(len=19) :: ch19, kbid, res19
    real(kind=8) :: r8b
    complex(kind=8) :: c16b
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    if (affeType .eq. 'CHAMP') then
!     -----------------------------
        ch19 = dsName
!
! ----- Create objects for global components (catalog) <=> local components (field)
!
        call cmpcha(ch19, cmp_name, cata_to_field, field_to_cata, nb_cmp)
        do k = 1, nb_cmp
            if (cmp_name(k) .ne. 'TEMP' .and. cmp_name(k) .ne. 'LAGR') then
                goto 7
            end if
        end do
        lautr = .false.
        goto 8
!
!
    else if (affeType .eq. 'EVOL') then
!     -----------------------------
        res19 = dsName
        call rsorac(res19, 'LONUTI', 0, r8b, kbid, &
                    c16b, r8b, kbid, nbordr, 1, &
                    ibid)
        call wkvect('&&AFVA01.NUME_ORDRE', 'V V I', nbordr(1), jordr)
        call rsorac(res19, 'TOUT_ORDRE', 0, r8b, kbid, &
                    c16b, r8b, kbid, zi(jordr), nbordr(1), &
                    ibid)
!
        do j = 1, nbordr(1)
            call rsexch('F', res19, fieldType, zi(jordr-1+j), ch19, &
                        iret)
            if (iret .eq. 0) then
!
! ------------- Create objects for global components (catalog) <=> local components (field)
!
                call cmpcha(ch19, cmp_name, cata_to_field, field_to_cata, nb_cmp)
                do k = 1, nb_cmp
                    if (cmp_name(k) .ne. 'TEMP' .and. cmp_name(k) .ne. 'LAGR') then
                        goto 7
                    end if
                end do
                AS_DEALLOCATE(vi=cata_to_field)
                AS_DEALLOCATE(vi=field_to_cata)
                AS_DEALLOCATE(vk8=cmp_name)
            end if
        end do
        call jedetr('&&AFVA01.NUME_ORDRE')
        lautr = .false.
        goto 8
!
!
    else
        write (6, *) affeType, dsName, fieldType
        ASSERT(.false.)
    end if
!
7   continue
    lautr = .true.
    goto 8
!
!
8   continue
!
    AS_DEALLOCATE(vi=cata_to_field)
    AS_DEALLOCATE(vi=field_to_cata)
    AS_DEALLOCATE(vk8=cmp_name)
!
    call jedema()
end subroutine
