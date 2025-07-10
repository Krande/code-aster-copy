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

subroutine copy_field_with_numbering(fieldin, fieldout, mesh, nume_equa, &
                                     base, typc, nequa)
    implicit none

#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/vtcrea.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtzero.h"
#include "asterfort/utmess.h"

    character(len=19), intent(in) :: fieldin
    character(len=19), intent(in) :: fieldout
    character(len=19), intent(in) :: mesh
    character(len=19), intent(in) :: nume_equa
    character(len=1), intent(in) :: base
    character(len=1), optional, intent(in) :: typc
    integer(kind=8), optional, intent(in) :: nequa
!_______________________________________________________________________
!     Copies field fieldin into a new field fieldout, converted to apply
!     numequa as numbering
!_______________________________________________________________________
!     Input :
!       fieldin   : name of input field (unchanged)
!       mesh      : name of mesh over which both fields are defined
!       nume_equa : name of target equation numbering
!       base      : into 'V', 'G' memory zone to write the copied field
!       typc (optional)  : real or complex type
!       nequa (optional) : total number of equations for target field
!     Input/output:
!       fieldout  : name of output field/ output field
!_______________________________________________________________________
!
    character(len=1) :: type_sca
    character(len=24) :: crefe(2)
    integer(kind=8):: iret, neq

    if (present(nequa)) then
        neq = nequa
    else
        call dismoi('NB_EQUA', nume_equa, 'NUME_EQUA', repi=neq)
    end if
    if (present(typc)) then
        type_sca = typc
    else
        ! For now retrieved from EquationNumbering
        ! To be seen if it should be retrieved from input fields
        call dismoi('TYPE_SCA', nume_equa, 'NUME_EQUA', repk=type_sca)
    end if

    crefe(1) = mesh
    crefe(2) = nume_equa
    call vtcrea(fieldout, crefe, base, type_sca, neq)
    if (type_sca == 'R') then
        call vtzero(fieldout)
    end if
    call vtcopy(fieldin, fieldout, iret)
!
    if (iret .ne. 0) then
        call utmess('A', 'UTILITAI_24')
    end if

end subroutine
