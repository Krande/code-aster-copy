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

subroutine dtmclean(sd_dtm_)
    implicit none

!
!
! person_in_charge: hassan.berro@edf.fr
!
! dtmclean : call all subroutines that need to be called to clean their
!            own data
!
!       sd_dtm_          : dtm data structure
!       sd_nl_           : nl  data structure
!
! =======================================================================

#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/dtmget.h"

!     1. Input / output arguments
    character(len=*), intent(in)  :: sd_dtm_
    character(len=8)                    :: sd_dtm, sd_nl

    integer(kind=8)                             :: nbnli

    call jemarq()

    sd_dtm = sd_dtm_

    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)

    if (nbnli .gt. 0) then

        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)

    end if

    call jedema()

end subroutine
