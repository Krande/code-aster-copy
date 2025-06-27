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
! ==================================================================================================
!
! Module for the contact (main algorithm)
!
! ==================================================================================================
!
module contactAlgo_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: addContactLigrel
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/gettco.h"
#include "asterfort/jeexin.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! addContactLigrel
!
! Add LIGREL from contact to list of LIGREL
!
! In  sdcontZ           : contact datastructure from DEFI_CONT/DEFI_CONTACT
! In  virtualCellZ      : virtual cells for contact solving
! IO  nbLigr            : number of LIGREL in list
! Ptr listLigr          : list of LIGREL
!
! --------------------------------------------------------------------------------------------------
    subroutine addContactLigrel(sdcontZ, virtualCellZ, nbLigr, listLigr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: sdcontZ, virtualCellZ
        integer(kind=8), intent(inout) :: nbLigr
        character(len=24), pointer :: listLigr(:)
! ----- Locals
        character(len=16) :: typeSD
        character(len=8) :: sdcont
        character(len=24) :: ligrelName, virtualCell
        aster_logical :: lDefiContact, lDefiCont, lParallelCont
        integer(kind=8) :: contForm, nbLigrCont, nbLigrNew, iret
        aster_logical :: lContCont, lContLac, lContDisc, lUnil
        character(len=24), pointer :: listLigrSave(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        sdcont = sdcontZ
        virtualCell = virtualCellZ

! ----- Detect original command
        call gettco(sdcont, typeSD)
        lDefiContact = ASTER_FALSE
        lDefiCont = ASTER_FALSE
        if (typeSD .eq. "CHAR_CONTACT") then
            lDefiContact = ASTER_TRUE
        elseif (typeSD .eq. "CHAR_CONT" .or. typeSD .eq. "CHAR_FROT" .or. &
                typeSD .eq. "LIGREL_CP") then
            lDefiCont = ASTER_TRUE
        else
            WRITE (6, *) "sdcont:", sdcont
            WRITE (6, *) "typeSD:", typeSD
            ASSERT(ASTER_FALSE)
        end if

        lParallelCont = typeSD .eq. "LIGREL_CP"

! ----- LIGREL
        nbLigrCont = 0
        lContCont = ASTER_FALSE
        lContLac = ASTER_FALSE
        lContDisc = ASTER_FALSE
        lUnil = ASTER_FALSE
        if (lDefiContact) then
! --------- Contact formulation
            contForm = cfdisi(sdcont(1:8)//'.CONTACT', 'FORMULATION')
            lContDisc = contForm .eq. 1
            lContCont = contForm .eq. 2
            lUnil = contForm .eq. 4
            lContLac = contForm .eq. 5

! --------- LIGREL from CONTINUE or LAC contact
            if (lContCont .or. lContLac) then
                ligrelName = sdcont(1:8)//'.CONT.LIGRE'
                call jeexin(ligrelName(1:19)//'.LIEL', iret)
                if (iret .gt. 0) then
                    nbLigrCont = nbLigrCont+1
                end if
                ligrelName = virtualCell
                call jeexin(ligrelName(1:19)//'.LIEL', iret)
                if (iret .gt. 0) then
                    nbLigrCont = nbLigrCont+1
                end if
            end if

        elseif (lDefiCont) then
            if (lParallelCont) then
                ligrelName = sdcont(1:8)
            else
                ligrelName = sdcont(1:8)//'.CONT.LIGRE'
            end if
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigrCont = nbLigrCont+1
            end if
            ligrelName = virtualCell
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigrCont = nbLigrCont+1
            end if

        elseif (lContDisc .or. lUnil) then
            nbLigrCont = 0

        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Update length of list
        nbLigrNew = nbLigr+nbLigrCont
        if (nbLigrNew .gt. nbLigr) then
            AS_ALLOCATE(vk24=listLigrSave, size=nbLigr)
            listLigrSave(1:nbLigr) = listLigr(1:nbLigr)
            AS_DEALLOCATE(vk24=listLigr)
            AS_ALLOCATE(vk24=listLigr, size=nbLigrNew)
            listLigr(1:nbLigr) = listLigrSave(1:nbLigr)
            AS_DEALLOCATE(vk24=listLigrSave)
        end if

! ----- Add virtual elements from CONTINUE or LAC contact
        if (lContCont .or. lContLac) then
            ligrelName = sdcont(1:8)//'.CONT.LIGRE'
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigr = nbLigr+1
                listLigr(nbLigr) = ligrelName
            end if
            ligrelName = virtualCell
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigr = nbLigr+1
                listLigr(nbLigr) = ligrelName
            end if
        elseif (lDefiCont) then
            if (lParallelCont) then
                ligrelName = sdcont(1:8)
            else
                ligrelName = sdcont(1:8)//'.CONT.LIGRE'
            end if
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigr = nbLigr+1
                listLigr(nbLigr) = ligrelName
            end if
            ligrelName = virtualCell
            call jeexin(ligrelName(1:19)//'.LIEL', iret)
            if (iret .gt. 0) then
                nbLigr = nbLigr+1
                listLigr(nbLigr) = ligrelName
            end if
        elseif (lContDisc .or. lUnil) then
            nbLigrCont = 0
        else
            ASSERT(ASTER_FALSE)
        end if

        ASSERT(nbLigr .eq. nbLigrNew)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module contactAlgo_module
