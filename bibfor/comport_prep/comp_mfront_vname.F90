! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
    integer, intent(in) :: model_dim
    integer, intent(in) :: nbVariMeca
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
    integer :: nbVariMFront, iVariType, iVari, iTens, variSize, leng
    character(len=16) :: vari_name, variName
    character(len=80), pointer :: variNameList(:) => null()
    integer, pointer :: variSizeList(:) => null()
    character(len=2), parameter :: cmpv_name(6) = (/'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'/)
    character(len=2), parameter :: cmpt_name(9) = (/ &
                                   'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8'/)
!
! --------------------------------------------------------------------------------------------------
!
    if (nbVariMeca .ne. 0) then
        call mgis_get_number_of_isvs(extern_addr, nbVariMFront)
        if (nbVariMFront .eq. 0) then
            iVari = 1
            infoVari(iVari) = 'VIDE'
        else
            AS_ALLOCATE(vk80=variNameList, size=nbVariMFront)
            AS_ALLOCATE(vi=variSizeList, size=nbVariMFront)
            call mgis_get_isvs(extern_addr, variNameList)
            call mgis_get_isvs_sizes(extern_addr, variSizeList)
            iVari = 0
            do iVariType = 1, nbVariMFront
                variName = variNameList(iVariType) (1:16)
                variSize = variSizeList(iVariType)
                leng = lxlgut(variName)
                ! scalar, vector, tensor
                ASSERT(variSize .eq. 1 .or. variSize .eq. 2*model_dim .or. variSize .eq. 9)
                if (variSize .eq. 1) then
                    infoVari(iVari+1) = variName
                elseif (variSize .eq. 2*model_dim) then
                    do iTens = 1, 2*model_dim
                        if (leng .le. 14) then
                            vari_name = variName(1:leng)//cmpv_name(iTens)
                        else
                            vari_name = variName(1:14)//cmpv_name(iTens)
                        end if
                        infoVari(iVari+iTens) = vari_name
                    end do
                elseif (variSize .eq. 9) then
                    do iTens = 1, 9
                        if (leng .le. 14) then
                            vari_name = variName(1:leng)//cmpt_name(iTens)
                        else
                            vari_name = variName(1:14)//cmpt_name(iTens)
                        end if
                        infoVari(iVari+iTens) = vari_name
                    end do
                end if
                iVari = iVari+variSize
            end do
            AS_DEALLOCATE(vk80=variNameList)
            AS_DEALLOCATE(vi=variSizeList)
        end if
        ASSERT(nbVariMeca .eq. iVari)
    end if
!
end subroutine
