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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romTableChck(tablReduCoor, lTablFromResu, nbModeIn, nbSnapIn_, nbStoreIn_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbexve.h"
#include "asterfort/tbGetListPara.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_TablReduCoor), intent(in) :: tablReduCoor
    aster_logical, intent(in) :: lTablFromResu
    integer(kind=8), intent(in) :: nbModeIn
    integer(kind=8), optional, intent(in) :: nbSnapIn_, nbStoreIn_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Check parameters of table for the reduced coordinates
!
! --------------------------------------------------------------------------------------------------
!
! In  tablReduCoor     : table for reduced coordinates
! In  lTablFromResu    : flag if table is in results datastructure
! In  nbModeIn         : number of modes in base
! In  nbSnapIn         : number of snapshots
! In  nbStoreIn        : number of storing index
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbPara, tablNbLine, nbMode, nbSnap, nbStore
    integer(kind=8) :: iPara, iLine
    character(len=24) :: tablName
    character(len=24), pointer :: paraName(:) => null()
    character(len=24), pointer :: paraType(:) => null()
    integer(kind=8), parameter :: nbParaRequired = 5
    character(len=24), parameter :: paraNameRequired(nbParaRequired) = (/ &
                                    'COOR_REDUIT', 'INST       ', &
                                    'NUME_MODE  ', 'NUME_ORDRE ', &
                                    'NUME_SNAP  '/)
    character(len=8), parameter :: paraTypeRequired(nbParaRequired) = (/ &
                                   'R', 'R', &
                                   'I', 'I', &
                                   'I'/)
    integer(kind=8), pointer :: snapNume(:) => null()
    integer(kind=8), pointer :: storeNume(:) => null()
    integer(kind=8) :: nbVale
    aster_logical :: lTablFromUser, lTablExist
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM15_1')
    end if
!
! - Get parameters
!
    lTablFromUser = tablReduCoor%lTablFromUser
    if (lTablFromUser) then
        tablName = tablReduCoor%tablUserName
    else
        tablName = tablReduCoor%tablResu%tablName
    end if
!
! - Is table exists ?
!
    lTablExist = lTablFromResu .or. lTablFromUser
    if (lTablExist) then
        if (lTablFromUser .and. lTablFromResu) then
            call utmess('F', 'ROM15_24')
        end if
    else
        call utmess('F', 'ROM15_23')
    end if
!
! - Get parameters in existing table
!
    call tbGetListPara(tablName, nbPara, paraName, paraType, tablNbLine)
!
! - Check number of parameters
!
    if (nbPara .ne. nbParaRequired) then
        call utmess('F', 'ROM15_27')
    end if
!
! - Check name/type of parameters
!
    do iPara = 1, nbPara
        if (paraName(iPara) .ne. paraNameRequired(iPara)) then
            call utmess('F', 'ROM15_27')
        end if
        if (paraType(iPara) .ne. paraTypeRequired(iPara)) then
            call utmess('F', 'ROM15_27')
        end if
    end do
    AS_DEALLOCATE(vk24=paraType)
    AS_DEALLOCATE(vk24=paraName)
!
! - Number of lines
!
    if (tablNbLine .eq. 0) then
        call utmess('F', 'ROM15_28')
    end if
!
! - Check number of snapshots
!
    if (present(nbSnapIn_)) then
        call tbexve(tablName, 'NUME_SNAP', '&&NUMESNAP', 'V', nbVale)
        call jeveuo('&&NUMESNAP', 'E', vi=snapNume)
        nbSnap = 0
        do iLine = 1, tablNbLine
            if (snapNume(iLine) .gt. nbSnap) then
                nbSnap = snapNume(iLine)
            end if
        end do
        if (nbSnap .ne. nbSnapIn_) then
            call utmess('F', 'ROM15_29')
        end if
        nbMode = tablNbLine/nbSnap
    end if
!
! - Check number of storing index
!
    if (present(nbStoreIn_)) then
        call tbexve(tablName, 'NUME_ORDRE', '&&NUMESTOREP', 'V', nbVale)
        call jeveuo('&&NUMESTOREP', 'E', vi=storeNume)
        nbStore = 0
        do iLine = 1, tablNbLine
            if (storeNume(iLine) .gt. nbStore) then
                nbStore = storeNume(iLine)
            end if
        end do
        if (nbStore .ne. nbStoreIn_) then
            call utmess('F', 'ROM15_30')
        end if
        nbMode = tablNbLine/nbStore
    end if
!
! - Check number of modes
!
    if (nbMode .ne. nbModeIn) then
        call utmess('F', 'ROM15_31')
    end if
!
end subroutine
