! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_mfront_vname(nbVariMeca , &
                             libr_name  , subr_name  , model_mfront, model_dim,&
                             infoVari)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/lxlgut.h"
#include "asterc/mfront_get_number_of_internal_state_variables.h"
#include "asterc/mfront_get_internal_state_variables.h"
#include "asterc/mfront_get_internal_state_variables_types.h"
!
integer, intent(in) :: nbVariMeca
character(len=255), intent(in) :: libr_name, subr_name
character(len=16), intent(in) :: model_mfront
integer, intent(in) :: model_dim
character(len=16), pointer :: infoVari(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Name of internal variables for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  nbVariMeca       : number of internal state variables for mechanical part of behaviour
! In  libr_name        : name of library
! In  subr_name        : name of comportement in library
! In  model_mfront     : type of modelisation MFront
! In  model_dim        : dimension of modelisation (2D or 3D)
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbVariType, iVariType, iVari, iTens, leng
    character(len=16) :: vari_name, variName, variType
    character(len=80), pointer :: variNameList(:) => null()
    character(len=80), pointer :: variTypeList(:) => null()
    character(len=2), parameter :: cmpv_name(6) = (/'XX','YY','ZZ','XY','XZ','YZ'/)
    character(len=2), parameter :: cmpt_name(9) = (/'F0','F1','F2','F3','F4','F5','F6','F7','F8'/)
!
! --------------------------------------------------------------------------------------------------
!
    call mfront_get_number_of_internal_state_variables(libr_name   , subr_name,&
                                                       model_mfront, nbVariType)
    if ( nbVariMeca .ne. 0 ) then
        AS_ALLOCATE(vk80 = variNameList, size = nbVariType)
        AS_ALLOCATE(vk80 = variTypeList, size = nbVariType)
        call mfront_get_internal_state_variables(libr_name, subr_name,&
                                                 model_mfront, variNameList,&
                                                 nbVariType)
        call mfront_get_internal_state_variables_types(libr_name, subr_name,&
                                                       model_mfront, variTypeList)
        iVari = 0
        do iVariType = 1, nbVariType
            variName = variNameList(iVariType)(1:16)
            variType = variTypeList(iVariType)(1:16)
            leng     = lxlgut(variName)
            if (variType .eq. 'scalar') then
                iVari = iVari + 1
                infoVari(iVari) = variName

            elseif (variType .eq. 'vector') then
                do iTens = 1, 2*model_dim
                    if (leng .le. 14) then
                        vari_name = variName(1:leng)//cmpv_name(iTens)
                    else
                        vari_name = variName(1:14)//cmpv_name(iTens)
                    endif
                    iVari = iVari + 1
                    infoVari(iVari) = vari_name
                end do

            elseif (variType .eq. 'tensor') then
                do iTens = 1, 9
                    if (leng .le. 14) then
                        vari_name = variName(1:leng)//cmpt_name(iTens)
                    else
                        vari_name = variName(1:14)//cmpt_name(iTens)
                    endif
                    iVari = iVari + 1
                    infoVari(iVari) = vari_name
                end do

            else
                ASSERT(ASTER_FALSE)

            endif
        end do
        AS_DEALLOCATE(vk80 = variNameList)
        AS_DEALLOCATE(vk80 = variTypeList)   
        ASSERT(nbVariMeca .eq. iVari)
    endif
!
end subroutine
