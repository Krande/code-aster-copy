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
subroutine nmdoch_wrap(listLoadZ, jvBaseZ, &
                       ligrel_slav0, ligrel_cont0)
!
    use listLoad_type
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmdoch.h"
!
    character(len=*), intent(in) :: listLoadZ
    character(len=*), intent(in) :: jvBaseZ
    character(len=*), intent(in) :: ligrel_slav0, ligrel_cont0
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics - Read parameters
!
! Get loads information and create datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad          : name of datastructure for list of loads
! In  jvBase            : JEVEUX base where to create objects
! In  ligrel_slav0      : late elements for slave side (contact)
! In  ligrel_cont0      : late elements for contact
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: jvBase
    character(len=24) :: listLoad
    integer :: nbLoadContact
    character(len=19) :: ligrel_slav, ligrel_cont
    type(ListLoad_Prep) :: listLoadPrep
    character(len=8), parameter :: funcCste = '&&NMDOME'
!
! --------------------------------------------------------------------------------------------------
!
    listLoad = listLoadZ
    ligrel_slav = ligrel_slav0
    ligrel_cont = ligrel_cont0
    jvBase = jvBaseZ
    nbLoadContact = 0
    if (ligrel_slav .ne. ' ') then
        nbLoadContact = nbLoadContact+1
    end if
    if (ligrel_cont .ne. ' ') then
        nbLoadContact = nbLoadContact+1
    end if

! - Standard loads
    call nmdoch(listLoadPrep, listLoad, jvBase)

! - "Loads" for contact
    if (nbLoadContact .ne. 0) then
        listLoadPrep%contactLigrelDefi = ligrel_slav
        listLoadPrep%contactLigrelSolv = ligrel_cont
        call addContactElements(listLoadPrep, &
                                jvBase, funcCste, &
                                listLoad)
    end if

! - Debug
    !call listLoadDebug(listLoad)
!
end subroutine
