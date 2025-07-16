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

subroutine pj3dcoSupInf(modelZ, meshNbNode, &
                        nbCellMast, cellMast, &
                        nbNodeSlav, nodeSlav, geomSlavJv, &
                        corre1, corre2, corre3, &
                        chnormZ, thickness)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pj3dco.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: modelZ
    integer(kind=8), intent(in):: meshNbNode
    integer(kind=8), intent(in) :: nbCellMast, cellMast(nbCellMast)
    integer(kind=8), intent(in) :: nbNodeSlav, nodeSlav(nbNodeSlav)
    character(len=24), intent(in) :: geomSlavJv
    character(len=16), intent(in) :: corre1, corre2, corre3
    character(len=*), intent(in) :: chnormZ
    real(kind=8), intent(in) :: thickness
!
! --------------------------------------------------------------------------------------------------
!
! BUT : CALCULER LES SD CORRE1 ET COORE2 UTILISEES POUR :
!       LIAISON_MAIL + TYPE_RACCORD='COQUE_MASSIF'
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model
    integer(kind=8) :: iNodeSlav, nodeNume, iDime, nbCmp
    character(len=19) :: chnorm
    character(len=19), parameter :: csnorm = '&&CALIR3.CSNORM'
    integer(kind=8), pointer :: cnsd(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    aster_logical, pointer :: cnsl(:) => null()
    real(kind=8), pointer :: geomSlav(:) => null()
    real(kind=8), pointer :: thickNormal(:) => null()
    aster_logical, parameter :: lDistMaxi = ASTER_FALSE
    real(kind=8), parameter :: distMaxi = 0.d0, distAlarm = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call jeveuo(geomSlavJv, 'E', vr=geomSlav)
    model = modelZ
    chnorm = chnormZ

! - Acces to normals
    call cnocns(chnorm, 'V', csnorm)
    call jeveuo(csnorm//'.CNSD', 'L', vi=cnsd)
    call jeveuo(csnorm//'.CNSL', 'L', vl=cnsl)
    call jeveuo(csnorm//'.CNSV', 'L', vr=cnsv)
    nbCmp = cnsd(2)
    ASSERT(nbCmp .eq. 3)
    ASSERT(nbCmp .eq. 3)

! - Compute normal vector with thickness
    call wkvect(corre3, 'V V R', 3*meshNbNode, vr=thickNormal)
    do iNodeSlav = 1, nbNodeSlav
        nodeNume = nodeSlav(iNodeSlav)
        do iDime = 1, 3
            if (.not. cnsl(3*(nodeNume-1)+iDime)) then
                call utmess('F', 'CHAMPS_2', sk=chnorm)
            end if
            thickNormal(3*(nodeNume-1)+iDime) = cnsv(3*(nodeNume-1)+iDime)*thickness
        end do
    end do

! - Compute projection for superior plane
    do iNodeSlav = 1, nbNodeSlav
        nodeNume = nodeSlav(iNodeSlav)
        do iDime = 1, 3
            geomSlav(3*(nodeNume-1)+iDime) = geomSlav(3*(nodeNume-1)+iDime)+ &
                                             thickNormal(3*(nodeNume-1)+iDime)/2.d0
        end do
    end do
    call pj3dco('PARTIE', model, model, nbCellMast, cellMast, &
                nbNodeSlav, nodeSlav, ' ', geomSlavJv, corre1, &
                lDistMaxi, distMaxi, distAlarm)

! - Compute projection for inferior plane
    do iNodeSlav = 1, nbNodeSlav
        nodeNume = nodeSlav(iNodeSlav)
        do iDime = 1, 3
            geomSlav(3*(nodeNume-1)+iDime) = geomSlav(3*(nodeNume-1)+iDime)- &
                                             thickNormal(3*(nodeNume-1)+iDime)
        end do
    end do
    call pj3dco('PARTIE', model, model, nbCellMast, cellMast, &
                nbNodeSlav, nodeSlav, ' ', geomSlavJv, corre2, &
                lDistMaxi, distMaxi, distAlarm)

! - Push back geometry
    do iNodeSlav = 1, nbNodeSlav
        nodeNume = nodeSlav(iNodeSlav)
        do iDime = 1, 3
            geomSlav(3*(nodeNume-1)+iDime) = geomSlav(3*(nodeNume-1)+iDime)+ &
                                             thickNormal(3*(nodeNume-1)+iDime)/2.d0
        end do
    end do
!
    call detrsd('CHAMP', csnorm)
    call jedema()
end subroutine
