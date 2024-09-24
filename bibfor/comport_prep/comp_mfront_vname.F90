! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine comp_mfront_vname(extern_addr, model_dim, nbVariMeca, infoVari)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/lxlgut.h"
#include "asterc/mgis_get_number_of_isvs.h"
#include "asterc/mgis_get_isvs.h"
#include "asterc/mgis_get_isvs_sizes.h"
!
    character(len=16), intent(in) :: extern_addr
    integer, intent(in) :: model_dim, nbVariMeca
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
! In  extern_addr      : MGIS address
! In  model_dim        : dimension of modelisation (2D or 3D)
! In  nbVariMeca       : number of internal state variables for mechanical part of behaviour
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbVariMGIS, iVariType, iVari, iCmp, variSizeMGIS, leng
    character(len=16) :: variName, variNameMGIS
    character(len=80), pointer :: variNameList(:) => null()
    integer, pointer :: variSizeList(:) => null()
    character(len=2), parameter :: cmpNameVoigt(6) = &
                                   (/'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'/)
    character(len=2), parameter :: cmpNameTensor(9) = &
                                   (/'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8'/)
    character(len=2), parameter :: cmpNameVector(3) = &
                                   (/'N ', 'T1', 'T2'/)
!
! --------------------------------------------------------------------------------------------------
!
    if (nbVariMeca .ne. 0) then
        call mgis_get_number_of_isvs(extern_addr, nbVariMGIS)
        if (nbVariMGIS .eq. 0) then
            iVari = 1
            infoVari(iVari) = 'VIDE'
        else
            AS_ALLOCATE(vk80=variNameList, size=nbVariMGIS)
            AS_ALLOCATE(vi=variSizeList, size=nbVariMGIS)
            call mgis_get_isvs(extern_addr, variNameList)
            call mgis_get_isvs_sizes(extern_addr, variSizeList)
            iVari = 0
            do iVariType = 1, nbVariMGIS
                variNameMGIS = variNameList(iVariType) (1:16)
                variSizeMGIS = variSizeList(iVariType)
                leng = lxlgut(variNameMGIS)
                if (variSizeMGIS .eq. 1) then
                    infoVari(iVari+1) = variNameMGIS
                elseif (variSizeMGIS .eq. 2*model_dim) then
                    do iCmp = 1, 2*model_dim
                        if (leng .le. 14) then
                            variName = variNameMGIS(1:leng)//cmpNameVoigt(iCmp)
                        else
                            variName = variNameMGIS(1:14)//cmpNameVoigt(iCmp)
                        end if
                        infoVari(iVari+iCmp) = variName
                    end do
                elseif (variSizeMGIS .eq. model_dim) then
                    do iCmp = 1, model_dim
                        if (leng .le. 14) then
                            variName = variNameMGIS(1:leng)//cmpNameVector(iCmp)
                        else
                            variName = variNameMGIS(1:14)//cmpNameVector(iCmp)
                        end if
                        infoVari(iVari+iCmp) = variName
                    end do
                elseif (variSizeMGIS .eq. 9) then
                    do iCmp = 1, 9
                        if (leng .le. 14) then
                            variName = variNameMGIS(1:leng)//cmpNameTensor(iCmp)
                        else
                            variName = variNameMGIS(1:14)//cmpNameTensor(iCmp)
                        end if
                        infoVari(iVari+iCmp) = variName
                    end do
                else
                    ASSERT(ASTER_FALSE)
                end if
                iVari = iVari+variSizeMGIS
            end do
            AS_DEALLOCATE(vk80=variNameList)
            AS_DEALLOCATE(vi=variSizeList)
        end if
        ASSERT(nbVariMeca .eq. iVari)
    end if
!
end subroutine
