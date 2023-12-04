! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine nmdini(keywf, list_inst, tole, &
                  nb_inst, l_init_noexist, nume_ini)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
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
    integer, intent(in) :: nb_inst
    aster_logical, intent(out) :: l_init_noexist
    integer, intent(out) :: nume_ini
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Index of initial time
!
! --------------------------------------------------------------------------------------------------
!
! In  keywf            : factor keyword
! In  list_inst        : list of times from INCREMENT/LIST_INST
! In  tole             : tolerance to search time
! In  nb_inst          : number of time steps in list
! Out l_init_noexist   : .true. if initial time doesn't exist in list of times
! Out nume_ini         : index of initial time
!
! --------------------------------------------------------------------------------------------------
!
    integer :: n1, n2
    real(kind=8) :: inst
    real(kind=8), pointer :: v_list_inst(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_ini = 0
    inst = -1.d0
    l_init_noexist = .false.
!
! - Acces to list of times
!
    call jeveuo(list_inst, 'L', vr=v_list_inst)
!
! - Get keywords
!
    call getvis(keywf, 'NUME_INST_INIT', iocc=1, scal=nume_ini, nbret=n1)
    call getvr8(keywf, 'INST_INIT', iocc=1, scal=inst, nbret=n2)
!
! - No NUME_INST_INIT/INST_INIT
!
    if (n1 .eq. 0) then
        if (n2 .eq. 0) then
            nume_ini = 0.0
        else
            call utacli(inst, v_list_inst, nb_inst, tole, nume_ini)
            if (nume_ini .lt. 0) then
                call utmess('F', 'DISCRETISATION_89', sr=inst)
            end if
        end if
    end if
!
! - Checks
!
    ASSERT(nume_ini .ge. 0)
    ASSERT(nume_ini .le. nb_inst)
!
end subroutine
