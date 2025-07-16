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

subroutine vrcomp_chck_rela(mesh, nbCell, &
                            comporCurr, comporPrev, &
                            ligrelCurr, ligrelPrev, &
                            comp_comb_1, comp_comb_2, verbose, &
                            newBehaviourOnCell, inconsistentBehaviour, &
                            l_modif_vari)
!
    use mesh_module, only: getGroupsFromCell
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cesexi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbCell
    character(len=19), intent(in) :: comporCurr, comporPrev
    character(len=19), intent(in) :: ligrelCurr, ligrelPrev
    character(len=48), intent(in) :: comp_comb_1, comp_comb_2
    aster_logical, intent(in) :: verbose
    aster_logical, intent(out) :: newBehaviourOnCell, inconsistentBehaviour, l_modif_vari

!
! --------------------------------------------------------------------------------------------------
!
! Check compatibility of comportments
!
! Check if comportments are the same (or compatible)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh          : name of mesh
! In  nbCell       : number of elements for current comportment
! In  comporCurr : reduced field for current comportment
! In  comporPrev : reduced field for previous comportment
! In  ligrelCurr   : current LIGREL
! In  ligrelPrev   : previous LIGREL
! In  comp_comb_1   : list of comportments can been mixed with each other
! In  comp_comb_2   : list of comportments can been mixed with all other ones
! Out newBehaviourOnCell    : .true. if not the same number of Gauss points
! Out inconsistentBehaviour  : .true. if not the same relation
! Out l_modif_vari  : .true. to change the structure of internal variables field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iCell
    aster_logical :: lCellCurr, lCellPrev
    integer(kind=8) :: iadp, iadm
    integer(kind=8) :: idxPrev, idxCurr
    character(len=16) :: relaCompPrev, relaCompCurr
    character(len=8) :: cellName
    integer(kind=8), pointer :: repePrev(:) => null()
    integer(kind=8), pointer :: repeCurr(:) => null()
    integer(kind=8) :: jvCeslPrev, jvCesdPrev
    character(len=16), pointer :: cesvPrev(:) => null()
    integer(kind=8) :: jvCeslCurr, jvCesdCurr
    character(len=24) :: valk(7)
    character(len=16), pointer :: cesvCurr(:) => null()
    character(len=24) :: groupCell(4)
    integer(kind=8) :: nbGroupCell
!
! --------------------------------------------------------------------------------------------------
!
    l_modif_vari = ASTER_FALSE
    newBehaviourOnCell = ASTER_FALSE
    inconsistentBehaviour = ASTER_FALSE

! - Access to LIGREL
    call jeveuo(ligrelCurr//'.REPE', 'L', vi=repeCurr)
    call jeveuo(ligrelPrev//'.REPE', 'L', vi=repePrev)

! - Acces to reduced CARTE on current comportement
    call jeveuo(comporCurr//'.CESD', 'L', jvCesdCurr)
    call jeveuo(comporCurr//'.CESV', 'L', vk16=cesvCurr)
    call jeveuo(comporCurr//'.CESL', 'L', jvCeslCurr)

! - Acces to reduced CARTE on previous comportement
    call jeveuo(comporPrev//'.CESD', 'L', jvCesdPrev)
    call jeveuo(comporPrev//'.CESV', 'L', vk16=cesvPrev)
    call jeveuo(comporPrev//'.CESL', 'L', jvCeslPrev)

! - Check on mesh
    do iCell = 1, nbCell
        lCellPrev = repePrev(2*(iCell-1)+1) .gt. 0
        lCellCurr = repeCurr(2*(iCell-1)+1) .gt. 0
        call cesexi('C', jvCesdPrev, jvCeslPrev, iCell, 1, 1, 1, iadm)
        call cesexi('C', jvCesdCurr, jvCeslCurr, iCell, 1, 1, 1, iadp)
        if (iadp .gt. 0) then
            relaCompCurr = cesvCurr(iadp)
            if (iadm .le. 0) then
                if (verbose) then
                    cellName = int_to_char8(iCell)
                    valk = ' '
                    valk(1) = cellName
                    call getGroupsFromCell(mesh, iCell, groupCell, nbGroupCell)
                    valk(1:nbGroupCell) = groupCell(1:nbGroupCell)
                    call utmess('I', 'COMPOR6_10', nk=6, valk=valk)
                end if
                newBehaviourOnCell = ASTER_TRUE
            else
                relaCompPrev = cesvPrev(iadm)

! ------------- Same comportement
                if (relaCompPrev .eq. relaCompCurr) then
                    goto 10
                else

! ----------------- These behaviours can been mixed ?
                    idxPrev = index(comp_comb_1, relaCompPrev)
                    idxCurr = index(comp_comb_1, relaCompCurr)
                    if ((idxPrev .gt. 0) .and. (idxCurr .gt. 0)) then
                        goto 10
                    end if

! ----------------- Comportements can been always mixed
                    idxPrev = index(comp_comb_2, relaCompPrev)
                    idxCurr = index(comp_comb_2, relaCompCurr)
                    if ((idxPrev .gt. 0) .or. (idxCurr .gt. 0)) then
                        l_modif_vari = ASTER_TRUE
                        goto 10
                    end if

! ----------------- Comportements cannot been mixed
                    if (lCellCurr .and. lCellPrev) then
                        if (verbose) then
                            cellName = int_to_char8(iCell)
                            valk = ' '
                            valk(1) = cellName
                            call getGroupsFromCell(mesh, iCell, groupCell, nbGroupCell)
                            valk(1:nbGroupCell) = groupCell(1:nbGroupCell)
                            call utmess('I', 'COMPOR6_10', nk=6, valk=valk)
                            valk = ' '
                            call getGroupsFromCell(mesh, iCell, groupCell, nbGroupCell)
                            call utmess('I', 'COMPOR6_11', nk=2, valk=[relaCompPrev, relaCompCurr])
                        end if
                        inconsistentBehaviour = ASTER_TRUE
                    end if
                end if
            end if
10          continue
        end if
    end do
!
end subroutine
