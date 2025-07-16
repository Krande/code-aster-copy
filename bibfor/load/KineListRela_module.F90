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
! Module to manage list of kinematic relations
!
! ==================================================================================================
!
module KineListRela_module
! ==================================================================================================
    use KineListRela_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: kineListRelaCreate, kineListRelaDelete, kineListRelaSave, kineListRelaInit
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/imprel.h"
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! kineListRelaCreate
!
! Create object for list of linear relations
!
! In  relaType         : type of relation
!                         'implicit': no second member (RHS)
!                         'explicit': with a second membre (RHS)
! In  nbTermMaxi       : maximum number of term for each linear relations
! In  listLineRela     : JEVEUX name of datastructure for list of linear relations
! Out kineListRela     : object for list of linear relations
! In  coefImpoType     : type of RHS (for explicit case)
!
! --------------------------------------------------------------------------------------------------
    subroutine kineListRelaCreate(relaTypeZ, nbTermMaxi, listLineRela, kineListRela, coefImpoType_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: relaTypeZ
        integer(kind=8), intent(in) :: nbTermMaxi
        character(len=19), intent(in) :: listLineRela
        type(KINE_LIST_RELA), intent(out) :: kineListRela
        character(len=4), intent(in), optional :: coefImpoType_
! - Local
        character(len=16) :: relaType
!   ------------------------------------------------------------------------------------------------
!
        relaType = relaTypeZ
        kineListRela%relaType = relaTypeZ
        kineListRela%listLineRela = listLineRela
        kineListRela%nbTermMaxi = nbTermMaxi

        if (relaType .eq. 'Implicit') then
! ----- Only real !
            kineListRela%coefImpoType = 'REEL'
            kineListRela%coefMultType = 'REEL'
            AS_ALLOCATE(vr=kineListRela%coefMultReal, size=nbTermMaxi)
            AS_ALLOCATE(vk8=kineListRela%nodeName, size=nbTermMaxi)
            AS_ALLOCATE(vk8=kineListRela%dofName, size=nbTermMaxi)
            AS_ALLOCATE(vr=kineListRela%LCSVale, size=nbTermMaxi*3)
            AS_ALLOCATE(vi=kineListRela%LCSType, size=nbTermMaxi)
        elseif (relaType .eq. 'Explicit') then
            if (coefImpoType_ .eq. 'REEL') then
                kineListRela%coefImpoType = 'REEL'
            elseif (coefImpoType_ .eq. 'FONC') then
                kineListRela%coefImpoType = 'FONC'
            else
                ASSERT(ASTER_FALSE)
            end if
            kineListRela%coefMultType = 'REEL'
            AS_ALLOCATE(vr=kineListRela%coefMultReal, size=nbTermMaxi)
            AS_ALLOCATE(vk8=kineListRela%nodeName, size=nbTermMaxi)
            AS_ALLOCATE(vk8=kineListRela%dofName, size=nbTermMaxi)
            AS_ALLOCATE(vr=kineListRela%LCSVale, size=nbTermMaxi*3)
            AS_ALLOCATE(vi=kineListRela%LCSType, size=nbTermMaxi)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineListRelaDelete
!
! Delete object for list of linear relations
!
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineListRelaDelete(kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
!   ------------------------------------------------------------------------------------------------
!
        AS_DEALLOCATE(vk8=kineListRela%nodeName)
        AS_DEALLOCATE(vk8=kineListRela%dofName)
        AS_DEALLOCATE(vr=kineListRela%coefMultReal)
        AS_DEALLOCATE(vr=kineListRela%LCSVale)
        AS_DEALLOCATE(vi=kineListRela%LCSType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineListRelaSave
!
! Save list of linear relations in object
!
! In  title            : title for debug
! In  nbTerm           : number of terms in linear relation
! In  kineListRela     : object for list of linear relations
! In  epsiDebg         : tolerance to print debug informations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineListRelaSave(titleZ, nbTerm, kineListRela, epsiDebg_)
! - Parameters
        character(len=*), intent(in) :: titleZ
        integer(kind=8), intent(in) :: nbTerm
        type(KINE_LIST_RELA), intent(in) :: kineListRela
        aster_logical, optional, intent(in) :: epsiDebg_
! - Local
        real(kind=8), parameter :: realZero = 0.d0
        complex(kind=8), parameter :: cplxZero = dcmplx(0.d0, 0.d0)
        character(len=8), parameter :: funcZero = ' '
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbTerm .le. kineListRela%nbTermMaxi)
        ASSERT(kineListRela%coefMultType .ne. 'COMP')
        ASSERT(kineListRela%coefImpoType .ne. 'COMP')

        if (kineListRela%relaType .eq. 'Implicit') then
! ----- Save list of relations
            call afrela(kineListRela%coefMultReal, [cplxZero], &
                        kineListRela%dofName, kineListRela%nodeName, &
                        kineListRela%LCSType, kineListRela%LCSVale, &
                        nbTerm, realZero, cplxZero, funcZero, &
                        kineListRela%coefMultType, kineListRela%coefImpoType, &
                        kineListRela%coefMultTole, kineListRela%listLineRela)

! ----- Print linear relation (for debug)
            if (present(epsiDebg_)) then
                call imprel(titleZ, nbTerm, &
                            kineListRela%coefMultReal, kineListRela%dofName, &
                            kineListRela%nodeName, &
                            realZero, kineListRela%coefMultTole)
            else
                call imprel(titleZ, nbTerm, &
                            kineListRela%coefMultReal, kineListRela%dofName, &
                            kineListRela%nodeName, &
                            realZero)
            end if

        elseif (kineListRela%relaType .eq. 'Explicit') then
! ----- Save list of relations
            call afrela(kineListRela%coefMultReal, [cplxZero], &
                        kineListRela%dofName, kineListRela%nodeName, &
                        kineListRela%LCSType, kineListRela%LCSVale, &
                        nbTerm, kineListRela%coefImpoReal, cplxZero, kineListRela%coefImpoFunc, &
                        kineListRela%coefMultType, kineListRela%coefImpoType, &
                        kineListRela%coefMultTole, kineListRela%listLineRela)

        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineListRelaInit
!
! Initialisation of list of relations
!
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineListRelaInit(kineListRela)
! - Parameters
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
!   ------------------------------------------------------------------------------------------------
!
        kineListRela%nodeName = ' '
        kineListRela%dofName = ' '
        kineListRela%coefMultReal = 0.d0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module KineListRela_module
