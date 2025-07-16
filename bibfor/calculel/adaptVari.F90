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
subroutine adaptVari(comporZ, variZ, ligrelZ)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cestas.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    character(len=*), intent(in) :: comporZ, variZ, ligrelZ
!
! --------------------------------------------------------------------------------------------------
!
! Modif VARI_ELGA field with new behaviours
!
! --------------------------------------------------------------------------------------------------
!
! In  compor       : map to describe behaviours
! In  vari         : internal variable field
! In  ligrel       : finite element descriptor
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: variS = '&&VRCOM2.CESV1'
    character(len=19), parameter :: newVariS = '&&VARCOM2.CESV2'
    character(len=19), parameter :: newVari = '&&VRCOM2.NEWVARI'
    character(len=19), parameter :: comporS = '&&VRCOM2.COTO'
    character(len=19), parameter :: comporSR = '&&VRCOM2.COPP'
    character(len=1), parameter :: base = 'G'
    integer(kind=8) :: iad1, iad2, nbCell, nbpg2, nbsp1, nbsp2, nbVari2, ipg, isp, iVari
    integer(kind=8) :: nbVari1
    integer(kind=8) :: iCell, iret
    integer(kind=8) :: iadp, jcoppl, jcoppd, jcoppv
    integer(kind=8) :: action
    integer(kind=8) :: jcev1d, jcev1l
    integer(kind=8) :: jcev2d, jcev2l, nncp, ibid
    real(kind=8), pointer :: cev1v(:) => null()
    real(kind=8), pointer :: cev2v(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Change VARI_ELGA field in simple field
    call celces(variZ, 'V', variS)
    call cestas(variS)
    call jeveuo(variS//'.CESD', 'L', jcev1d)
    call jeveuo(variS//'.CESV', 'L', vr=cev1v)
    call jeveuo(variS//'.CESL', 'L', jcev1l)
    nbCell = zi(jcev1d-1+1)

! - Create field of right dimensions
    call alchml(ligrelZ, 'RAPH_MECA', 'PVARIPR', 'V', newVari, iret, comporZ)
    call celces(newVari, 'V', newVariS)
    call detrsd('CHAM_ELEM', newVari)
    call jeveuo(newVariS//'.CESD', 'L', jcev2d)
    call jeveuo(newVariS//'.CESV', 'E', vr=cev2v)
    call jeveuo(newVariS//'.CESL', 'E', jcev2l)

! - Extract field of RELA_NAME if behaviours' map
    call carces(comporZ, 'ELEM', ' ', 'V', comporS, 'A', iret)
    call cesred(comporS, 0, [0], 1, 'RELCOM', 'V', comporSR)
    call detrsd('CHAM_ELEM_S', comporS)
    call jeveuo(comporSR//'.CESD', 'L', jcoppd)
    call jeveuo(comporSR//'.CESV', 'L', jcoppv)
    call jeveuo(comporSR//'.CESL', 'L', jcoppl)

! - Loop on cells
    do iCell = 1, nbCell
! ----- Sizes of new field
        nbpg2 = zi(jcev2d-1+5+4*(iCell-1)+1)
        nbsp2 = zi(jcev2d-1+5+4*(iCell-1)+2)
        nbVari2 = zi(jcev2d-1+5+4*(iCell-1)+3)

! ----- This cell doesn't exist in new LIGREL
        if (nbsp2 .eq. 0) cycle
        call cesexi('C', jcoppd, jcoppl, iCell, 1, 1, 1, iadp)
        if (iadp .le. 0) cycle

! ----- Sizes of previous field
        nbsp1 = zi(jcev1d-1+5+4*(iCell-1)+2)
        nbVari1 = zi(jcev1d-1+5+4*(iCell-1)+3)
!
!       -- PARFOIS LE COMPORTEMENT EST AFFECTE SUR LES MAILLES
!          DE BORD ALORS QUE CES ELEMENTS N'ONT PAS DE VARIABLES
!          INTERNES (I.E. ILS IGNORENT RAPH_MECA).
!          ON NE VEUT PAS FAIRE D'ERREUR <F> :
        if ((nbsp1 .eq. 0) .and. (nbVari1 .eq. 0)) cycle
        ASSERT(nbsp2 .eq. nbsp1)
!
        if (nbVari1 .eq. nbVari2) then
            action = 1
        else
            action = 2
        end if
!
        do ipg = 1, nbpg2
            do isp = 1, nbsp2
                do iVari = 1, nbVari2
                    call cesexi('S', jcev2d, jcev2l, iCell, ipg, isp, iVari, iad2)
                    ASSERT(iad2 .gt. 0)
                    zl(jcev2l-1+iad2) = ASTER_TRUE
                    if (action .eq. 1 .or. &
                        ((action .eq. 2) .and. iVari .le. nbVari1)) then
                        call cesexi('S', jcev1d, jcev1l, iCell, ipg, isp, iVari, iad1)
                        ASSERT(iad1 .gt. 0)
                        cev2v(iad2) = cev1v(iad1)
                    else
                        cev2v(iad2) = 0.d0
                    end if
                end do
            end do
        end do
    end do

! - Create new field (after projection)
    call detrsd('CHAM_ELEM', variZ)
    call cescel(newVariS, ligrelZ, 'RAPH_MECA', 'PVARIMR', 'OUI', nncp, base, variZ, 'F', ibid)

! - Clean
    call detrsd('CHAM_ELEM_S', variS)
    call detrsd('CHAM_ELEM_S', newVariS)
    call detrsd('CHAM_ELEM_S', comporSR)

    call jedema()
end subroutine
