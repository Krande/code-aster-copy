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
subroutine nmdoct(listLoad, ds_contact)
!
    use NonLin_Datastructure_type
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "LoadTypes_type.h"
!
    character(len=19), intent(in) :: listLoad
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Initializations
!
! Prepare list of loads (and late elements) for contact
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad         : list of loads
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbLoadIden = 1
    character(len=24) :: listLoadIden(LOAD_NBIDEN_MAXI)
    integer(kind=8) :: iLoad, iret, nbLoad
    character(len=8) :: ligrel_link, lag12
    character(len=24) :: loadName24, loadNameJv
    character(len=8) :: loadFunc, funcCste
    character(len=1), parameter :: jvBase = "V"
    character(len=4), parameter :: phenom = "MECA"
    character(len=24), pointer :: listLoadName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    if (ds_contact%l_contact) then

! ----- Check if single Lagrange multiplier is not used
        loadNameJv = listLoad(1:19)//'.LCHA'
        call jeveuo(loadNameJv, 'L', vk24=listLoadName)
        call getNbLoadsFromList(listLoad, nbLoad)
        do iLoad = 1, nbLoad
            loadName24 = listLoadName(iLoad)
            call jeexin(loadName24(1:19)//'.LGRF', iret)
            if (iret .ne. 0) then
                call dismoi('TYPE_LAGR', loadName24(1:19), 'LIGREL', repk=lag12)
                if (lag12 .eq. 'LAG1') then
                    call utmess('F', 'MECANONLINE_5')
                end if
            end if
        end do

! ----- Prepare constant function
        funcCste = '&&NMDOCT'
        call createUnitFunc(funcCste, 'V', loadFunc)

! ----- Add list of linear relation
        if (ds_contact%l_dof_rela) then
            ligrel_link = ds_contact%ligrel_dof_rela
            listLoadIden(1) = "DIRI_CSTE"
            call addLoadToList(phenom, listLoad, jvBase, &
                               ligrel_link, loadFunc, &
                               nbLoadIden, listLoadIden)
        end if
    end if

    !call listLoadDebug(listLoad)
!
end subroutine
