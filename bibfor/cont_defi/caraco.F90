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

subroutine caraco(sdcont, keywf, cont_form, nb_cont_zone)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/caralv.h"
#include "asterfort/cazoco.h"
#include "asterfort/cazocp.h"
#include "asterfort/cazofm.h"
!
    character(len=8), intent(in) :: sdcont
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: cont_form
    integer(kind=8), intent(in) :: nb_cont_zone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get parameters of contact
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  keywf            : factor keyword to read
! In  cont_form        : formulation of contact
! In  nb_cont_zone     : number of zones of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_zone
!
! --------------------------------------------------------------------------------------------------
!

! - Get method of contact
    call cazofm(sdcont, keywf, cont_form, nb_cont_zone)

! - Get parameters (not depending on contact zones)
    call cazocp(sdcont)

! - Get parameters (depending on contact zones)
    do i_zone = 1, nb_cont_zone
        call cazoco(sdcont, keywf, cont_form, i_zone, &
                    nb_cont_zone)
    end do

! - Set automatic parameters
    call caralv(sdcont, nb_cont_zone, cont_form)
!
end subroutine
