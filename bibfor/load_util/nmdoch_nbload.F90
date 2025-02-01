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
subroutine nmdoch_nbload(l_zero_allowed, nbLoad, &
                         factorKeyword)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/getvid.h"
#include "asterfort/utmess.h"
!
    aster_logical, intent(in) :: l_zero_allowed
    integer, intent(out) :: nbLoad
    character(len=16), intent(out) :: factorKeyword
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics/Thermics - Read parameters
!
! Get number of loads for loads datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  l_zero_allowed  : .true. if we can create "zero-load" list of loads datastructure
! Out nbLoad          : number of loads for list of loads
! Out factorKeyword   : factor keyword to read loads
!                      'None' => no load to read
!                      'EXCIT' or ' ' => depending on command
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iFactorKeyword, iret_cable, iret_cable_cine, nocc, nbFactorKeyword
    character(len=8) :: load_name
!
! --------------------------------------------------------------------------------------------------
!
    nbLoad = 0
    factorKeyword = 'None'

! - Detect factor keyword for loads
    if (getexm('EXCIT', 'CHARGE') .eq. 1) then
        factorKeyword = 'EXCIT'
    end if
    if (getexm(' ', 'CHARGE') .eq. 1) then
        factorKeyword = ' '
    end if

! - Count number of loads
    if (factorKeyword .eq. 'None') then
        nbLoad = 0
    elseif (factorKeyword .eq. '  ') then
        call getvid(factorKeyword, 'CHARGE', iocc=0, nbret=nocc)
        nbLoad = abs(nocc)
    elseif (factorKeyword .eq. 'EXCIT') then
        call getfac(factorKeyword, nbFactorKeyword)
        do iFactorKeyword = 1, nbFactorKeyword
            call getvid(factorKeyword, 'CHARGE', iocc=iFactorKeyword, &
                        scal=load_name, nbret=nocc)
!
! ------------- For DEFI_CABLE_BP: count load only if kinematic
! ------------- (because Neumann is not load but initial stress)
!
            if (nocc .eq. 1) then
                call jeexin(load_name//'.CHME.SIGIN.VALE', iret_cable)
                if (iret_cable .eq. 0) then
                    nbLoad = nbLoad+1
                else
                    call jeexin(load_name//'.CHME.CIMPO.DESC', iret_cable_cine)
                    if (iret_cable_cine .ne. 0) then
                        nbLoad = nbLoad+1
                    end if
                end if
            end if
        end do
    else
        ASSERT(.false.)
    end if

! - No loads is allowed ?
    if (nbLoad .eq. 0) then
        if (.not. l_zero_allowed) then
            call utmess('F', 'CHARGES_2')
        end if
    end if
!
end subroutine
