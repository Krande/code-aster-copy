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
subroutine vrcomp_chck_cmp(mesh, nbCell, &
                           comporCurrZ, &
                           comporCurr, comporPrev, &
                           variRedu, comp_comb_2, &
                           ligrelCurr, ligrelPrev, &
                           verbose, &
                           nbSpgDifferent, nbVariDifferent, l_modif_vari)
!
    use mesh_module, only: getGroupsFromCell
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbCell
    character(len=*), intent(in) :: comporCurrZ
    character(len=19), intent(in) :: comporCurr, comporPrev
    character(len=19), intent(in) :: variRedu
    character(len=48), intent(in) :: comp_comb_2
    character(len=19), intent(in) :: ligrelCurr, ligrelPrev
    aster_logical, intent(in) :: verbose
    aster_logical, intent(out) :: nbSpgDifferent, nbVariDifferent, l_modif_vari
!
! --------------------------------------------------------------------------------------------------
!
! Check compatibility of comportments
!
! Check if elements have the same number of internal variables and Gauss-subpoints
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh          : name of mesh
! In  nbCell        : number of elements for current comportment
! In  comporCurr    : current comportment
! In  comporCurr    : reduced field for current comportment
! In  comporPrev    : reduced field for previous comportment
! In  variRedu      : reduced field for internal variable
! In  comp_comb_2   : list of comportments can been mixed with all other ones
! In  ligrelCurr    : current LIGREL
! In  ligrelPrev    : previous LIGREL
! Out nbSpgDifferent   : .true. if not the same number of Gauss-subpoints
! Out nbVariDifferent   : .true. if not the same number of components
! Out l_modif_vari  : .true. to change the structure of internal variables field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iad1, iad2, iadp, iadm
    integer(kind=8) :: iCell, k, vali(3)
    aster_logical :: lCellCurr, lCellPrev
    integer(kind=8) :: idxPrev, idxCurr
    aster_logical :: all_is_zero
    integer(kind=8) :: nbPgPrev, nbSpgPrev, nbVariPrev
    integer(kind=8) :: nbSpgCurr, nbVariCurr
    character(len=16) :: relaCompPrev, relaCompCurr
    character(len=8) :: cellName
    character(len=19) :: dcel
    integer(kind=8), pointer :: repePrev(:) => null()
    integer(kind=8), pointer :: repeCurr(:) => null()
    integer(kind=8) :: jvDcelCesd, jvDcelCesl
    integer(kind=8), pointer :: dcelCesv(:) => null()
    integer(kind=8) :: jvCompCurrCesl, jvCompCurrCesd
    character(len=16), pointer :: compCurrCesv(:) => null()
    integer(kind=8) :: jvCompPrevCesd, jvCompPrevCesl
    character(len=16), pointer :: compPrevCesv(:) => null()
    integer(kind=8) :: jvVariCesd, jvVariCesl
    real(kind=8), pointer :: variCesv(:) => null()
    character(len=24) :: groupCell(4), valk(7)
    integer(kind=8) :: nbGroupCell
!
! --------------------------------------------------------------------------------------------------
!
    l_modif_vari = ASTER_FALSE
    nbVariDifferent = ASTER_FALSE
    nbSpgDifferent = ASTER_FALSE

! - Access to LIGREL
    call jeveuo(ligrelCurr//'.REPE', 'L', vi=repeCurr)
    call jeveuo(ligrelPrev//'.REPE', 'L', vi=repePrev)

! - Acces to reduced field on internal variables
    call jeveuo(variRedu//'.CESD', 'L', jvVariCesd)
    call jeveuo(variRedu//'.CESV', 'L', vr=variCesv)
    call jeveuo(variRedu//'.CESL', 'L', jvVariCesl)

! - Access to reduced CARTE DCEL_I (see CESVAR) on current behaviour
    dcel = comporCurrZ
    call jeveuo(dcel//'.CESD', 'L', jvDcelCesd)
    call jeveuo(dcel//'.CESV', 'L', vi=dcelCesv)
    call jeveuo(dcel//'.CESL', 'L', jvDcelCesl)

! - Acces to reduced CARTE on current behaviour
    call jeveuo(comporCurr//'.CESD', 'L', jvCompCurrCesd)
    call jeveuo(comporCurr//'.CESV', 'L', vk16=compCurrCesv)
    call jeveuo(comporCurr//'.CESL', 'L', jvCompCurrCesl)

! - Acces to reduced CARTE on previous behaviour
    if (comporPrev .ne. ' ') then
        call jeveuo(comporPrev//'.CESD', 'L', jvCompPrevCesd)
        call jeveuo(comporPrev//'.CESV', 'L', vk16=compPrevCesv)
        call jeveuo(comporPrev//'.CESL', 'L', jvCompPrevCesl)
    end if

! - Check on mesh
    do iCell = 1, nbCell
        lCellPrev = repePrev(2*(iCell-1)+1) .gt. 0
        lCellCurr = repeCurr(2*(iCell-1)+1) .gt. 0

! ----- Access to number of "sub-points"
        call cesexi('C', jvDcelCesd, jvDcelCesl, iCell, 1, 1, 1, iad1)
! ----- Access to number of internal state variables
        call cesexi('C', jvDcelCesd, jvDcelCesl, iCell, 1, 1, 2, iad2)
! ----- Access to current behaviour
        call cesexi('C', jvCompCurrCesd, jvCompCurrCesl, iCell, 1, 1, 1, iadp)

! ----- No behaviour on this element -> next element
        if (iad1 .le. 0) then
            cycle
        end if

! ----- Number of Gauss points/components
        ASSERT(iad2 .gt. 0)
        nbSpgCurr = dcelCesv(iad1)
        nbVariCurr = dcelCesv(iad2)
        nbPgPrev = zi(jvVariCesd-1+5+4*(iCell-1)+1)
        nbSpgPrev = zi(jvVariCesd-1+5+4*(iCell-1)+2)
        nbVariPrev = zi(jvVariCesd-1+5+4*(iCell-1)+3)

! ----- Check number of Gauss sub-points
        if (nbSpgCurr .ne. 0 .and. nbSpgPrev .ne. 0) then
            if (nbSpgCurr .ne. nbSpgPrev) then
                if (verbose) then
                    cellName = int_to_char8(iCell)
                    valk = ' '
                    valk(1) = cellName
                    call getGroupsFromCell(mesh, iCell, groupCell, nbGroupCell)
                    valk(1:nbGroupCell) = groupCell(1:nbGroupCell)
                    call utmess('I', 'COMPOR6_10', nk=6, valk=valk)
                    vali(1) = nbSpgPrev
                    vali(2) = nbSpgCurr
                    call utmess('I', 'COMPOR6_12', sk=cellName, ni=2, vali=vali)
                end if
                nbSpgDifferent = ASTER_TRUE
            end if
        end if

! ----- Check number of internal state variables
        if (nbVariCurr .ne. nbVariPrev) then

! --------- No vari => skin element
            if ((nbVariPrev .eq. 0) .or. (nbVariCurr .eq. 0)) then
                cycle
            end if

! --------- Current comportement can been mixed -> no problem
            ASSERT(iadp .gt. 0)
            relaCompCurr = compCurrCesv(iadp)
            idxCurr = index(comp_comb_2, relaCompCurr)
            if (idxCurr .gt. 0) then
                l_modif_vari = ASTER_TRUE
                cycle
            end if

! --------- Previous comportement can been mixed -> no problem
            if (comporPrev .ne. ' ') then
                call cesexi('C', jvCompPrevCesd, jvCompPrevCesl, iCell, 1, 1, 1, iadm)
                ASSERT(iadm .gt. 0)
                relaCompPrev = compPrevCesv(iadm)
                idxPrev = index(comp_comb_2, relaCompPrev)
                if (idxPrev .gt. 0) then
                    l_modif_vari = ASTER_TRUE
                    cycle
                end if
            else
                if (nbVariPrev .eq. 1) then
                    call cesexi('C', jvVariCesd, jvVariCesl, iCell, 1, 1, 1, iad2)
                    ASSERT(iad2 .gt. 0)
                    all_is_zero = ASTER_TRUE
                    do k = 1, nbPgPrev*nbSpgCurr
                        if (variCesv(iad2+k-1) .ne. 0.d0) then
                            all_is_zero = ASTER_FALSE
                        end if
                    end do
                    if (all_is_zero) then
                        l_modif_vari = ASTER_TRUE
                        cycle
                    end if
                end if
            end if
        end if

! ----- Not the same number of internal state variables
        if (nbVariCurr .ne. nbVariPrev) then
            l_modif_vari = ASTER_TRUE
            nbVariDifferent = ASTER_TRUE
            if (verbose) then
                cellName = int_to_char8(iCell)
                valk = ' '
                valk(1) = cellName
                call getGroupsFromCell(mesh, iCell, groupCell, nbGroupCell)
                valk(1:nbGroupCell) = groupCell(1:nbGroupCell)
                call utmess('I', 'COMPOR6_10', nk=6, valk=valk)
                vali(1) = nbVariPrev
                vali(2) = nbVariCurr
                call utmess('I', 'COMPOR6_13', sk=cellName, ni=2, vali=vali)
            end if
        end if
    end do
!
end subroutine
