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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine checkListOfGrpMa(mesh, listGrpMa, nbGrpMa, l_stop_local)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/existGrpMa.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in)    :: mesh
    character(len=*), intent(in)    :: listGrpMa(*)
    integer(kind=8), intent(in)             :: nbGrpMa
    aster_logical, intent(in)       :: l_stop_local
!
!---------------------------------------------------------------------------------------------------
!   But :
!     Check list of groups of cells
!
!   IN:
!     mesh       : name of the mesh
!     l_stop     : emmit Fatal error if a group of cells is not in the local mesh
!     listGrpMa  : list of groups of cells to check
!     nbGrpMa    : number of groups of cells in the list
!
!---------------------------------------------------------------------------------------------------
    character(len=24) :: GrpMaName, valk(2)
    integer(kind=8) :: iGrpMa
    aster_logical ::  l_exi_in_grp, l_exi_in_grp_p
!-----------------------------------------------------------------------
!
    call jemarq()
!
    do iGrpMa = 1, nbGrpMa
        GrpMaName = listGrpMa(iGrpMa)
        call existGrpMa(mesh(1:8), GrpMaName, l_exi_in_grp, l_exi_in_grp_p)
!
! --- The group of cells exists in the global mesh
!
        if (l_exi_in_grp_p) then
!
! --- The group of cells does not exist in the local mesh
!
            if (.not. l_exi_in_grp .and. l_stop_local) then
                valk(1) = GrpMaName
                valk(2) = mesh
                call utmess('F', 'MODELISA7_77', nk=2, valk=valk)
            end if
        else
            valk(1) = GrpMaName
            valk(2) = mesh
            call utmess('F', 'MODELISA7_77', nk=2, valk=valk)
        end if
    end do
!
    call jedema()
!
end subroutine
