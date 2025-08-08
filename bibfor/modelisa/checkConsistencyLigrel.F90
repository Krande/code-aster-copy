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

subroutine checkConsistencyLigrel(model, ligrel, answer)
    implicit none
!
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
!
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: ligrel
    aster_logical, intent(out) :: answer
! BUT :
! Check that ligrel is included in model
!     -----------------------------------------------------------------
    integer(kind=8) :: nbCell, iCell
    character(len=8) :: ma1, ma2
    character(len=19) :: moLigrel
    integer(kind=8), pointer :: typfemo(:) => null()
    integer(kind=8), pointer :: typfe(:) => null()
!----------------------------------------------------------------------

    answer = ASTER_TRUE
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=moLigrel)
    call dismoi('NOM_MAILLA', moLigrel, 'LIGREL', repk=ma1)
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma2)
!
    if (ma1 .ne. ma2) then
        answer = ASTER_FALSE
        go to 999
    end if
!
    call dismoi('NB_MA_MAILLA', ma1, 'MAILLAGE', repi=nbCell)
!
    call jeveuo(moLigrel//'.TYFE', 'L', vi=typfemo)
    call jeveuo(ligrel//'.TYFE', 'L', vi=typfe)
!
    do iCell = 1, nbCell
        if (typfe(iCell) .ne. 0) then
            if (typfe(iCell) .ne. typfemo(iCell)) then
                answer = ASTER_FALSE
                go to 999

            end if
        end if
    end do
!
999 continue
!
end subroutine
