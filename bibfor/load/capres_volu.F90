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
subroutine capres_volu(load, mesh, valeType, nbOccPresRep)
!
    use SolidShell_Mesh_module, only: setValueOnFace
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/getelem.h"
#include "asterfort/dismoi.h"
!
    character(len=4), intent(in) :: valeType
    character(len=8), intent(in) :: load, mesh
    integer(kind=8), intent(in) :: nbOccPresRep
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads PRES_REP on volumic elements (solid shells)
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : name of load
! In  mesh             : name of mesh
! In  valeType         : affected value type (real, complex or function)
! In  nbOccPresRep     : number of occurrences for PRES_REP
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywFact = 'PRES_REP'
    integer(kind=8), parameter :: nbCmp = 2
    character(len=8), parameter :: cmpName(nbCmp) = (/'PINF', 'PSUP'/)
    character(len=8), pointer :: mapCmpName(:) => null()
    real(kind=8), pointer :: mapCmpValeR(:) => null()
    character(len=8), pointer :: mapCmpValeK(:) => null()
    real(kind=8) :: presUserR
    character(len=8) :: presUserK
    integer(kind=8), pointer :: meshTypmail(:) => null()
    character(len=19) :: carte
    character(len=24), parameter :: cellSkinJv = '&&LIST_ELEM'
    integer(kind=8), pointer :: cellSkin(:) => null()
    integer(kind=8), pointer :: presCell(:) => null()
    real(kind=8), pointer :: presFaceR(:) => null()
    character(len=8), pointer :: presFaceK(:) => null()
    integer(kind=8) :: nbCellSkin, nbCrack, nbCellHexa9, nbCellMesh
    integer(kind=8) :: iocc, iCell, iCellSkin, cellVoluNume
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Name of CARTE
    carte = load//'.CHME.PRESS'
!
! - Pre-allocation of CARTE
!
!   NOTHING => already done in capres_skin
!
! - Set name of components
    call jeveuo(carte//'.NCMP', 'E', vk8=mapCmpName)
    mapCmpName(1) = cmpName(1)
    mapCmpName(2) = cmpName(2)

! - Access to values
    if (valeType .eq. 'REEL') then
        call jeveuo(carte//'.VALV', 'E', vr=mapCmpValeR)
    else
        call jeveuo(carte//'.VALV', 'E', vk8=mapCmpValeK)
    end if

! - Access to mesh
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', nbCellMesh)
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)

! - Number of HEXA9 cells
    nbCellHexa9 = 0
    do iCell = 1, nbCellMesh
        if (meshTypmail(iCell) .eq. MT_HEXA9) then
            nbCellHexa9 = nbCellHexa9+1
        end if
    end do
    ASSERT(nbCellHexa9 .gt. 0)

! - Prepare vector
    AS_ALLOCATE(vi=presCell, size=nbCellHexa9)
    AS_ALLOCATE(vr=presFaceR, size=2*nbCellHexa9)
    AS_ALLOCATE(vk8=presFaceK, size=2*nbCellHexa9)
    presFaceK = '&FOZERO'

! - Set values in CARTE
    do iocc = 1, nbOccPresRep

! ----- Get value from user
        presUserR = 0.d0
        presUserK = '&FOZERO'
        if (valeType .eq. 'REEL') then
            call getvr8(keywFact, 'PRES', iocc=iocc, scal=presUserR)
        else
            call getvid(keywFact, 'PRES', iocc=iocc, scal=presUserK)
        end if

! ----- No cracks
        call getvid(keywFact, 'FISSURE', iocc=iocc, nbval=0, nbret=nbCrack)
        if (nbCrack .ne. 0) then
            call utmess('F', 'CHARGES_7')
        end if

! ----- Get list of faces from user
        call getelem(mesh, keywFact, iocc, 'A', cellSkinJv, nbCellSkin)
        call jeveuo(cellSkinJv, 'L', vi=cellSkin)

! ----- Get volume elements from list of faces
        call setValueOnFace(mesh, &
                            presUserR, presUserK, &
                            nbCellSkin, cellSkin, &
                            presCell, presFaceR, presFaceK)
!
        call jedetr(cellSkinJv)
    end do

! - Set values
    do iCellSkin = 1, nbCellSkin
        cellVoluNume = presCell(iCellSkin)
        if (valeType .eq. 'REEL') then
            mapCmpValeR(1) = presFaceR(2*(iCellSkin-1)+1)
            mapCmpValeR(2) = presFaceR(2*(iCellSkin-1)+2)
        else
            mapCmpValeK(1) = presFaceK(2*(iCellSkin-1)+1)
            mapCmpValeK(2) = presFaceK(2*(iCellSkin-1)+2)
        end if
        if (cellVoluNume .ne. 0) then
            call nocart(carte, 3, nbCmp, mode='NUM', nma=1, &
                        limanu=[cellVoluNume])
        end if
    end do

! - Clean
    AS_DEALLOCATE(vi=presCell)
    AS_DEALLOCATE(vr=presFaceR)
    AS_DEALLOCATE(vk8=presFaceK)
!
    call jedema()
end subroutine
