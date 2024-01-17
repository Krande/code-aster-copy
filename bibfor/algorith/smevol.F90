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

subroutine smevol(resultName, nbStore, listStore, &
                  model, materialField, materialCoding, comporMeta, &
                  metaInitUser, numeFieldInit)
!
    use MetallurgyOperator_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: resultName
    integer, intent(in) :: nbStore
    integer, pointer :: listStore(:)
    character(len=8), intent(in) :: model, materialField
    character(len=24), intent(in) :: materialCoding
    character(len=24), intent(in) :: comporMeta
    integer, intent(in) :: numeFieldInit
    character(len=24), intent(in) :: metaInitUser
!
!   ------------------------------------------------------------------------------------------------
!
!     Compute META_ELNO
!
!   ------------------------------------------------------------------------------------------------
!
    integer :: iret, jvPara, iStore, numphi
    integer :: numeStore_0, numeStore_1, numeStore_2
    real(kind=8) :: time_0, time_1, time_2
    character(len=19) :: resultName19
    character(len=24) :: modelLigrel, resultField
    character(len=24) :: chftrc, metaIn, metaOut
    character(len=24) :: temp_0, temp_1, temp_2

!   ------------------------------------------------------------------------------------------------
!
    call jemarq()
    resultName19 = resultName

    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Prepare field to manage TRC curves in elementary computation
    call metaPrepTRCWorkingField(model, materialField, chftrc)

! - Prepare dynamic field from behaviour
    call detrsd('CHAM_ELEM_S', comporMeta)
    call cesvar(' ', comporMeta, modelLigrel, comporMeta)

! - Compute initial field for metallurgy if required and store it
    if (numeFieldInit .eq. 0) then
        numphi = 1
        numeStore_1 = listStore(1)
        call metaPrepareInitialState(resultName, numeStore_1, &
                                     model, materialCoding, comporMeta, metaInitUser)
        numeStore_2 = listStore(2)
        call metaPrepareInitialState(resultName, numeStore_2, &
                                     model, materialCoding, comporMeta, metaInitUser)
    else
        numphi = 0
        do iStore = 2, nbStore
            if (listStore(iStore) .eq. numeFieldInit) then
                numphi = iStore-1
            end if
        end do
    end if

! - Main loop to compute
    do iStore = 1, nbStore-2
!
        numeStore_0 = listStore(iStore)
        if (numeStore_0 .lt. listStore(numphi)) cycle

! ----- Compute metallurgy from time 1 to time 2
        numeStore_1 = listStore(iStore+1)
        numeStore_2 = listStore(iStore+2)

! ----- Get fields
        call rsexch('F', resultName, 'TEMP', numeStore_0, temp_0, iret)
        call rsexch('F', resultName, 'TEMP', numeStore_1, temp_1, iret)
        call rsexch('F', resultName, 'META_ELNO', numeStore_1, metaIn, iret)
        call rsexch('F', resultName, 'TEMP', numeStore_2, temp_2, iret)

! ----- Get times
        call rsadpa(resultName, 'L', 1, 'INST', numeStore_0, 0, sjv=jvPara)
        time_0 = zr(jvPara)
        call rsadpa(resultName, 'L', 1, 'INST', numeStore_1, 0, sjv=jvPara)
        time_1 = zr(jvPara)
        call rsadpa(resultName, 'L', 1, 'INST', numeStore_2, 0, sjv=jvPara)
        time_2 = zr(jvPara)

! ----- Compute option META_ELNO
        call metaCompMetaElno(model, materialCoding, &
                              comporMeta, chftrc, &
                              time_0, time_1, time_2, &
                              temp_0, temp_1, temp_2, &
                              metaIn, metaOut)

! ----- Save metallurgy field
        call rsexch(' ', resultName, 'META_ELNO', numeStore_2, resultField, iret)
        call copisd('CHAMP_GD', 'G', metaOut, resultField)
        call rsnoch(resultName, 'META_ELNO', numeStore_2)
        call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStore_2, sr=time_2)

! ----- Save behaviour map
        call rsexch(' ', resultName, 'COMPORMETA', numeStore_2, resultField, iret)
        call copisd('CHAMP_GD', 'G', comporMeta, resultField)
        call rsnoch(resultName, 'COMPORMETA', numeStore_2)
        call utmess('I', 'ARCHIVAGE_6', sk='COMPORMETA', si=numeStore_2, sr=time_2)
    end do

! - Cleaning
    call jedetr('&&SMEVOL_FTRC')
    call jedetr('&&SMEVOL_TRC')
    call detrsd('CARTE', chftrc)
!
    call jedema()
end subroutine
