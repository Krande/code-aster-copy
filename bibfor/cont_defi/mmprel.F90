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

subroutine mmprel(sdcont, mesh, slavElemLigr)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ajellt.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: sdcont, mesh
    character(len=19), intent(in) :: slavElemLigr
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Continue method - Create slave elements
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  model            : name of model
! In  mesh             : name of mesh
! In  slavElemLigr     : LIGREL for virtual elements (slave side)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lFricZone, l_veri
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=16) :: modeli
    character(len=16), parameter :: phenom = "MECANIQUE"
    integer(kind=8) :: jdecme, iContZone
    integer(kind=8) :: nbContZone, model_ndim, ntElemSlav
    integer(kind=8) :: iCellSlav, cellSlavNume, nbCellSlav
    aster_logical :: l_verif_all
    character(len=24), parameter :: listCellJv = '&&MMPREL.LISTE_MAILLES'
    integer(kind=8), pointer :: listCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Datastructure for contact definition
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)

! - Parameters
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    ntElemSlav = cfdisi(sdcont_defi, 'NTMAEC')
    nbContZone = cfdisi(sdcont_defi, 'NZOCO')
    l_verif_all = cfdisl(sdcont_defi, 'ALL_VERIF')

! - Add elements
    if (.not. l_verif_all) then
! ----- Create list of slave elements
        call wkvect(listCellJv, 'V V I', ntElemSlav, vi=listCell)

! ----- Set list of slave elements
        do iContZone = 1, nbContZone
! --------- Type of model
            lFricZone = mminfl(sdcont_defi, 'FROTTEMENT_ZONE', iContZone)
            l_veri = mminfl(sdcont_defi, 'VERIF', iContZone)
            if (model_ndim .eq. 2) then
                if (lFricZone) then
                    modeli = 'FRIC_SL_2D'
                else
                    modeli = 'CONT_SL_2D'
                end if
            else if (model_ndim .eq. 3) then
                if (lFricZone) then
                    modeli = 'FRIC_SL_3D'
                else
                    modeli = 'CONT_SL_3D'
                end if
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Type of model
            if (.not. l_veri) then
                nbCellSlav = mminfi(sdcont_defi, 'NBMAE', iContZone)
                jdecme = mminfi(sdcont_defi, 'JDECME', iContZone)
                ASSERT(nbCellSlav .le. ntElemSlav)
                do iCellSlav = 1, nbCellSlav
                    cellSlavNume = v_sdcont_mailco(jdecme+iCellSlav)
                    listCell(iCellSlav) = cellSlavNume
                end do
                call ajellt(slavElemLigr, mesh, nbCellSlav, listCell, &
                            phenom, modeli)
            end if
        end do
        call jedetr(listCellJv)
    end if
!
end subroutine
