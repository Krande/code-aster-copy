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

subroutine ltnotb(result, table_iden, table_name, iret_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/gnoms2.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
!
!
    character(len=*), intent(in) :: result
    character(len=*), intent(in) :: table_iden
    character(len=*), intent(out) :: table_name
    integer(kind=8), optional, intent(out) :: iret_
!
! --------------------------------------------------------------------------------------------------
!
! Results datastructure
!
! Get/Create object for table in results datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of results datastructure
! In  table_iden       : identifier of table
! Out table_name       : name of JEVEUX object for this table
! Out iret             : return code
!                         iret = 0 => this table exists
!                         iret = 1 => this table doesn't exist and we cannot create it
!                         iret = 2 => this table doesn't exist and we can create it
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: noojb
    integer(kind=8) :: nb_table_maxi, nb_table_curr, i_table, iret_obj, iret
    character(len=16) :: table_iden_k16
    character(len=16), pointer :: v_tabl_ltnt(:) => null()
    character(len=24), pointer :: v_tabl_ltns(:) => null()
    character(len=19) :: result_k19
    character(len=24) :: valk(2)
    aster_logical :: stop_error
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    iret = 2
    table_name = ' '
    result_k19 = result
    table_iden_k16 = table_iden
    stop_error = .true._1
    if (present(iret_)) then
        stop_error = .false._1
    end if
!
! - Get object for list of tables in datastructure
!
    call jeexin(result_k19//'.LTNT', iret_obj)
    if (iret_obj .eq. 0) then
        iret = 1
        if (stop_error) then
            valk(1) = table_iden
            valk(2) = result
            call utmess('F', 'TABLE0_37', nk=2, valk=valk)
        end if
        goto 999
    end if
!
! - Access to objects
!
    call jelira(result_k19//'.LTNT', 'LONMAX', nb_table_maxi)
    call jelira(result_k19//'.LTNT', 'LONUTI', nb_table_curr)
    call jeveuo(result_k19//'.LTNT', 'L', vk16=v_tabl_ltnt)
!
! - Look for this table in list
!
    do i_table = 1, nb_table_curr
        if (v_tabl_ltnt(i_table) .eq. table_iden_k16) then
            call jeveuo(result_k19//'.LTNS', 'L', vk24=v_tabl_ltns)
            table_name = v_tabl_ltns(i_table)
            iret = 0
            goto 999
        end if
    end do
!
! - Create this table in list
!
    iret = 2
    nb_table_curr = nb_table_curr+1
    if (nb_table_curr .gt. nb_table_maxi) then
        call juveca(result_k19//'.LTNT', nb_table_curr+6)
        call juveca(result_k19//'.LTNS', nb_table_curr+6)
    end if
    call jeecra(result_k19//'.LTNT', 'LONUTI', nb_table_curr)
    call jeecra(result_k19//'.LTNS', 'LONUTI', nb_table_curr)
!
    call jeveuo(result_k19//'.LTNT', 'E', vk16=v_tabl_ltnt)
    v_tabl_ltnt(nb_table_curr) = table_iden_k16
    noojb = result_k19(1:8)//'.TB000000  .TBBA'
    call gnoms2(noojb, 12, 17)
    call jeveuo(result_k19//'.LTNS', 'E', vk24=v_tabl_ltns)
    v_tabl_ltns(nb_table_curr) = noojb(1:17)
    table_name = v_tabl_ltns(nb_table_curr)
!
999 continue
!
    if (present(iret_)) then
        iret_ = iret
    end if
!
    call jedema()
end subroutine
