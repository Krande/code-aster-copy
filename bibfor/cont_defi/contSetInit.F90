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
subroutine contSetInit(sdcont, mesh, nbContZone)
!
    use mesh_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnbsf.h"
#include "asterfort/cfzone.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfr.h"
!
    character(len=8), intent(in) :: sdcont, mesh
    integer(kind=8), intent(in) :: nbContZone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Continue method - Set value for CONTACT_INIT='INTERPENETRE'
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  mesh             : name of mesh
! In  nbContZone       : number of zones of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iContZone, iContSurf, contInit, jdecma, jdecsl, iCell, zcmcf
    real(kind=8) :: contInitDist, edgeMin, edgeMax
    integer(kind=8) :: nbCell, nbCellMast, nbCellSlav
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco, sdcont_caracf
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    real(kind=8), pointer :: v_sdcont_caracf(:) => null()
    integer(kind=8), pointer :: listCellNume(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
    sdcont_caracf = sdcont_defi(1:16)//'.CARACF'
    zcmcf = cfmmvd('ZCMCF')

    do iContZone = 1, nbContZone
        contInit = mminfi(sdcont_defi, 'CONTACT_INIT', iContZone)
        contInitDist = mminfr(sdcont_defi, 'CONTACT_INIT_DIST', iContZone)
        call cfzone(sdcont_defi, iContZone, 'ESCL', iContSurf)
        call cfnbsf(sdcont_defi, iContSurf, 'MAIL', nbCellSlav, jdecsl)
        call cfzone(sdcont_defi, iContZone, 'MAIT', iContSurf)
        call cfnbsf(sdcont_defi, iContSurf, 'MAIL', nbCellMast, jdecma)
        nbCell = nbCellMast+nbCellSlav
        AS_ALLOCATE(vi=listCellNume, size=nbCell)
        do iCell = 1, nbCellSlav
            listCellNume(iCell) = v_sdcont_mailco(jdecsl+iCell)
        end do
        do iCell = 1, nbCellMast
            listCellNume(nbCellSlav+iCell) = v_sdcont_mailco(jdecma+iCell)
        end do
        if (contInit .eq. 2) then
            if (contInitDist .lt. 0.d0) then
                call compMinMaxEdges(mesh, edgeMin, edgeMax, &
                                     nbCell, listCellNume)
                call jeveuo(sdcont_caracf, 'E', vr=v_sdcont_caracf)
                v_sdcont_caracf(zcmcf*(iContZone-1)+17) = edgeMin
            end if
        end if
        AS_DEALLOCATE(vi=listCellNume)
    end do
!
end subroutine
