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
subroutine resuSaveParameters(resultName, resultType, &
                              model, caraElem, materField, listLoadResu, &
                              empiNumePlan_, empiSnapNb_, empiFieldType_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
    character(len=8), intent(in) :: resultName
    character(len=16), intent(in) :: resultType
    character(len=8), intent(in) :: model, caraElem, materField
    character(len=24), intent(in) :: listLoadResu
    integer(kind=8), optional, intent(in) :: empiNumePlan_, empiSnapNb_
    character(len=24), optional, intent(in) :: empiFieldType_
!
! --------------------------------------------------------------------------------------------------
!
! LIRE_RESU and CREA_RESU
!
! Save standard parameters in results datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  resultName       : name of results datastructure
! In  resultType       : type of results datastructure (EVOL_NOLI, EVOL_THER, )
! In  model            : name of model
! In  caraElem         : name of elementary characteristics field
! In  materField       : name of material field
! In  listLoadResu     : name of datastructure for loads
! In  empiNumeplan     : index of plane for empiric modes
! In  empiSnapNb       : number of snapshots for empiric modes
! In  empiFieldType    : type of field for empiric modes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: storeNb, jvPara, iStore, storeNume
    integer(kind=8), pointer :: storeList(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Get number of storing slots in results datastructure
    call rs_get_liststore(resultName, storeNb)
    if (storeNb .le. 0) then
        call utmess('F', 'RESULT2_97')
    end if

! - Get list of storing slots in results datastructure
    AS_ALLOCATE(vi=storeList, size=storeNb)
    call rs_get_liststore(resultName, storeNb, storeList)

! - Save material field
    if (materField .ne. ' ') then
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsadpa(resultName, 'E', 1, 'CHAMPMAT', storeNume, 0, sjv=jvPara)
            zk8(jvPara) = materField
        end do
    end if

! - Save elementary characteristics field
    if (caraElem .ne. ' ') then
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsadpa(resultName, 'E', 1, 'CARAELEM', storeNume, 0, sjv=jvPara)
            zk8(jvPara) = caraElem
        end do
    end if

! - Save model
    if (model .ne. ' ') then
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsadpa(resultName, 'E', 1, 'MODELE', storeNume, 0, sjv=jvPara)
            zk8(jvPara) = model
        end do
    end if

! - Save list of loads
    if (listLoadResu .ne. ' ') then
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsadpa(resultName, 'E', 1, 'EXCIT', storeNume, 0, sjv=jvPara)
            zk24(jvPara) = listLoadResu
        end do
    end if

! - Save parameters for MODE_EMPI
    if (resultType .eq. 'MODE_EMPI') then
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsadpa(resultName, 'E', 1, 'NUME_PLAN', storeNume, 0, sjv=jvPara)
            zi(jvPara) = empiNumePlan_
            call rsadpa(resultName, 'E', 1, 'NB_SNAP', storeNume, 0, sjv=jvPara)
            zi(jvPara) = empiSnapNb_
            call rsadpa(resultName, 'E', 1, 'NOM_CHAM', storeNume, 0, sjv=jvPara)
            zk24(jvPara) = empiFieldType_
        end do
    end if

! - Clean
    AS_DEALLOCATE(vi=storeList)
!
end subroutine
