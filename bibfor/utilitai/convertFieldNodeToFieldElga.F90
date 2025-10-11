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
subroutine convertFieldNodeToFieldElga(model, fieldNode, fieldElga)
!
    implicit none
!
#include "asterfort/alchml.h"
#include "asterfort/chpchd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/nopar2.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model
    character(len=*), intent(in) :: fieldNode, fieldElga
!
! --------------------------------------------------------------------------------------------------
!
! Utility
!
! Convert nodal field to elga field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'TOU_INI_ELGA'
    character(len=19), parameter :: fieldElemRefe = '&&ELGA.CELMOD'
    character(len=19) :: ligrel
    character(len=8) :: paraName, physName
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------

    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
!
! - Create new field on cell
!
    call dismoi('NOM_GD', fieldNode, 'CHAMP', repk=physName, arret='C', ier=iret)
    call nopar2(option, physName, 'OUT', paraName)
    call alchml(ligrel, option, paraName, 'V', fieldElemRefe, iret, ' ')
    if (iret .ne. 0) then
        call utmess('F', 'UTILITAI3_23', nk=3, valk=[ligrel, paraName, option])
    end if
!
! - Change support of field
!
    call chpchd(fieldNode, 'ELGA', fieldElemRefe, 'OUI', 'G', fieldElga, model)
    call detrsd('CHAMP', fieldElemRefe)
!
end subroutine
