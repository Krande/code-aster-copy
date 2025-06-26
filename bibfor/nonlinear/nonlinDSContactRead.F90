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
subroutine nonlinDSContactRead(ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/getvid.h"
!
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Contact management
!
! Read parameters for contact management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nocc
    character(len=8) :: sdcont
    character(len=16), parameter :: factorKeyword = 'CONTACT'
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_10')
    end if

! - Get name of datastructure from DEFI_CONTACT
    call getvid(' ', factorKeyword, scal=sdcont, nbret=nocc)
    ds_contact%l_contact = nocc .gt. 0
    if (nocc .ne. 0) then
        ds_contact%sdcont = sdcont
    end if
!
end subroutine
