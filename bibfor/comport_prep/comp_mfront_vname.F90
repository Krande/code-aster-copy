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
subroutine comp_mfront_vname(extern_addr, nbVariMeca, infoVari)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/codent.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterc/mgis_get_number_of_isvs.h"
#include "asterc/mgis_get_isvs.h"
#include "asterc/mgis_get_isvs_sizes.h"
#include "asterc/mgis_get_isvs_types.h"
!
    character(len=16), intent(in) :: extern_addr
    integer(kind=8), intent(in) :: nbVariMeca
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
! In  nbVariMeca       : number of internal state variables for mechanical part of behaviour
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
  integer(kind=8) :: nbVariMGIS, iVariType, iVari, iCmp, variSizeMGIS, leng, variTypeMGIS, lenTronca
    integer(kind=8) :: lenMaxi
    character(len=16) :: variName, variNameMGIS
    character(len=80), pointer :: variNameList(:) => null()
    integer(kind=8), pointer :: variSizeList(:) => null()
    integer(kind=8), pointer :: variTypeList(:) => null()
    character(len=2) :: nomk2
    character(len=2), parameter :: cmpNameSTensor(6) = &
                                   (/'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'/)
    character(len=2), parameter :: cmpNameTensor(9) = &
                                   (/'XX', 'YY', 'ZZ', 'XY', 'YX', 'XZ', 'ZX', 'YZ', 'ZY'/)
    character(len=2), parameter :: cmpNameVector(3) = &
                                   (/'X', 'Y', 'Z'/)
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
            AS_ALLOCATE(vi=variTypeList, size=nbVariMGIS)
            call mgis_get_isvs(extern_addr, variNameList)
            call mgis_get_isvs_sizes(extern_addr, variSizeList)
            call mgis_get_isvs_types(extern_addr, variTypeList)
            iVari = 0
            do iVariType = 1, nbVariMGIS
                variNameMGIS = variNameList(iVariType) (1:16)
                variSizeMGIS = variSizeList(iVariType)
                variTypeMGIS = variTypeList(iVariType)

                leng = lxlgut(variNameMGIS)
                lenTronca = min(leng, 14)

                if (variTypeMGIS .eq. MGIS_BV_SCALAR) then
                    lenMaxi = 16
                    infoVari(iVari+1) = variNameMGIS

                else if (variTypeMGIS .eq. MGIS_BV_VECTOR_1D .or. &
                         variTypeMGIS .eq. MGIS_BV_VECTOR_2D .or. &
                         variTypeMGIS .eq. MGIS_BV_VECTOR_3D .or. &
                         variTypeMGIS .eq. MGIS_BV_VECTOR) then
                    lenMaxi = 14
                    if (variSizeMGIS .le. 3) then
                        do iCmp = 1, variSizeMGIS
                            variName = variNameMGIS(1:lenTronca)//cmpNameVector(iCmp)
                            infoVari(iVari+iCmp) = variName
                        end do
                    else
                        call utmess('F', "MGIS1_20")
                    end if

                else if (variTypeMGIS .eq. MGIS_BV_STENSOR_1D .or. &
                         variTypeMGIS .eq. MGIS_BV_STENSOR_2D .or. &
                         variTypeMGIS .eq. MGIS_BV_STENSOR_3D .or. &
                         variTypeMGIS .eq. MGIS_BV_STENSOR) then
                    lenMaxi = 14
                    if (variSizeMGIS .le. 6) then
                        do iCmp = 1, variSizeMGIS
                            variName = variNameMGIS(1:lenTronca)//cmpNameSTensor(iCmp)
                            infoVari(iVari+iCmp) = variName
                        end do
                    else
                        call utmess('F', "MGIS1_21")
                    end if

                else if (variTypeMGIS .eq. MGIS_BV_TENSOR_1D .or. &
                         variTypeMGIS .eq. MGIS_BV_TENSOR_2D .or. &
                         variTypeMGIS .eq. MGIS_BV_TENSOR_3D .or. &
                         variTypeMGIS .eq. MGIS_BV_TENSOR) then
                    lenMaxi = 14
                    if (variSizeMGIS .le. 9) then
                        do iCmp = 1, variSizeMGIS
                            variName = variNameMGIS(1:lenTronca)//cmpNameTensor(iCmp)
                            infoVari(iVari+iCmp) = variName
                        end do
                    else
                        call utmess('F', "MGIS1_22")
                    end if

                else if (variTypeMGIS .eq. MGIS_BV_HIGHER_ORDER_TENSOR .or. &
                         variTypeMGIS .eq. MGIS_BV_ARRAY) then
                    lenMaxi = 14
                    if (variSizeMGIS .le. 99) then
                        do iCmp = 1, variSizeMGIS
                            call codent(iCmp, 'D0', nomk2)
                            variName = variNameMGIS(1:lenTronca)//nomk2
                            infoVari(iVari+iCmp) = variName
                        end do
                    else
                        call utmess('F', "MGIS1_23")
                    end if

                else
                    call utmess('F', "MGIS1_24")
                end if
                if (leng .gt. lenMaxi) then
                    call utmess('A', "MGIS1_19", nk=2, &
                                valk=[variNameMGIS, variNameMGIS(1:lenTronca)])
                end if
                iVari = iVari+variSizeMGIS
            end do
            AS_DEALLOCATE(vk80=variNameList)
            AS_DEALLOCATE(vi=variSizeList)
            AS_DEALLOCATE(vi=variTypeList)
        end if
        ASSERT(nbVariMeca .eq. iVari)
    end if
!
end subroutine
