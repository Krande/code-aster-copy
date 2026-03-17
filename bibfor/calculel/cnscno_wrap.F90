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

subroutine cnscno_wrap(cnsz, nume_equaz, prol0, basez, cnoz, kstop, iret)
!
    implicit none
!
#include "asterfort/cnscno.h"
#include "asterfort/global_numbering_compute.h"
#include "asterfort/global_numbering_communicate.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
    !
    character(len=*) :: cnsz, cnoz, basez, nume_equaz, prol0
    character(len=1) :: kstop
    integer(kind=8) :: iret
    aster_logical :: l_exi_nume
    character(len=19) :: nume_equa
!
    l_exi_nume = .false.
    if (nume_equaz .ne. ' ') then
        call exisd("NUME_EQUA", nume_equaz, iret)
        if (iret .eq. 1) then
            l_exi_nume = .true.
        end if
    end if
    call cnscno(cnsz, nume_equaz, prol0, basez, cnoz, kstop, iret, lprofconst=ASTER_FALSE)
    if (.not. l_exi_nume) then
        call dismoi('NUME_EQUA', cnoz, 'CHAM_NO', repk=nume_equa)
        ! create numbering
        call global_numbering_compute(nume_equa)
        ! communicate numbering
        call global_numbering_communicate(nume_equa)
    end if
end subroutine
