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

subroutine nmdifi(keywf, list_inst, tole, nb_inst, nume_end)
!
    implicit none
!
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utacli.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=16), intent(in) :: keywf
    character(len=19), intent(in) :: list_inst
    real(kind=8), intent(in) :: tole
    integer(kind=8), intent(in) :: nb_inst
    integer(kind=8), intent(out) :: nume_end
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Index of final time
!
! --------------------------------------------------------------------------------------------------
!
! In  keywf            : factor keyword
! In  list_inst        : list of times from INCREMENT/LIST_INST
! In  tole             : tolerance to search time
! In  nb_inst          : number of time steps in list
! Out nume_end         : index of final time
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n1, n2
    real(kind=8) :: inst
    real(kind=8), pointer :: v_list_inst(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_end = 0
!
! - Acces to list of times
!
    call jeveuo(list_inst, 'L', vr=v_list_inst)
!
! - Get keywords
!
    call getvis(keywf, 'NUME_INST_FIN', iocc=1, scal=nume_end, nbret=n1)
    call getvr8(keywf, 'INST_FIN', iocc=1, scal=inst, nbret=n2)
!
! - No NUME_INST_FIN/INST_FIN
!
    if (n1+n2 .eq. 0) then
        nume_end = nb_inst-1
    else if (n1 .eq. 0) then
        call utacli(inst, v_list_inst, nb_inst, tole, nume_end)
    end if

! - Checks
    if (nume_end .lt. 0 .or. nume_end .gt. (nb_inst-1)) then
        call utmess('F', 'DISCRETISATION_94')
    end if
!
end subroutine
