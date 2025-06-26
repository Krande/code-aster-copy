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

subroutine apntos(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apcaln.h"
#include "asterfort/apforc.h"
#include "asterfort/apvepa.h"
#include "asterfort/infdbg.h"
#include "asterfort/mmbouc.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Pairing - Node to segment
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, err_appa
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<Pairing> Node-to-segment pairing'
    end if
    err_appa = 0
!
! - Compute tangents
!
    call apcaln(mesh, ds_contact, err_appa)
!
! - Pairing by "brute" force
!
    call apforc(mesh, ds_contact, err_appa)

    if (err_appa .eq. 1) then
        call mmbouc(ds_contact, 'Geom', 'Set_Error')
    else
        call mmbouc(ds_contact, 'Geom', 'Set_NoError')
    end if
!
! - Check pairing
!
    call apvepa(ds_contact)
!
end subroutine
