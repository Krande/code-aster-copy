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
subroutine ccbcop(resultIn, resultOut, &
                  listStoreJv, nbStore, &
                  listOptionJv, nbOption)
!
    use postComp_module
    use postComp_type
    use result_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/ccfnrn.h"
#include "asterfort/ccfnrnLegacy.h"
#include "asterfort/deprecated_option.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rs_get_caraelem.h"
#include "asterfort/rs_get_model.h"
#include "asterfort/rscrsd.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbStore, nbOption
    character(len=8), intent(in) :: resultOut, resultIn
    character(len=19), intent(in) :: listStoreJv, listOptionJv
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! BOUCLE SUR LA LISTE D'OPTION ET APPEL A CALCOP
!
! --------------------------------------------------------------------------------------------------
!
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! In  nbStore           : number of storing indexes
! In  listStoreJv       : JEVEUX name object for list of storing indexes
! In  nbOption          : number of options to compute
! In  listOptionJv      : JEVEUX name object for list of options to compute
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, numeStore0, iOption, codret
    integer(kind=8), pointer :: listStore(:) => null()
    character(len=16), pointer :: listOption(:) => null()
    character(len=8) :: modelRefe, caraElemRefe, answer
    character(len=16) :: option, resultType
    aster_logical :: newResult, lTimeDist, lRDM
    type(POST_COMP) :: postComp
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Management of result datastructures
    call gettco(resultIn, resultType)
    call jeexin(resultOut//'           .DESC', iret)
    newResult = iret .eq. 0
    if ((resultIn .ne. resultOut) .and. (.not. newResult)) then
        call utmess('F', 'CALCCHAMP1_18')
    end if

! - New output result
    if (newResult) then
        call rscrsd('G', resultOut, resultType, nbStore)
        call titre()
    end if

! - Access to list of storing indexes
    call jeveuo(listStoreJv, 'L', vi=listStore)

! - Get first storing index
    if (nbStore .lt. 1) then
        call utmess("F", "CALCCHAMP1_1")
    end if
    numeStore0 = listStore(1)

! - Get reference model
    call rs_get_model(resultIn, numeStore0, modelRefe, codret)
    if (codret .lt. 0) then
        call utmess('F', 'CALCCHAMP1_44')
    end if

! - Detect structural elements
    call dismoi('EXI_RDM', modelRefe, 'MODELE', repk=answer)
    lRDM = answer(1:3) .eq. 'OUI'

! - Get reference CARA_ELEM
    call rs_get_caraelem(resultIn, numeStore0, caraElemRefe, codret)
    if (lRDM) then
        if (codret .lt. 0) then
            call utmess('A', 'CALCCHAMP1_45')
        end if
    end if

! - Copy parameters
    if (newResult) then
        call rsCopyPara(resultIn, resultOut, nbStore, listStore)
    end if

! - Access to list of options to compute
    call jeexin(listOptionJv, iret)
    if (iret .ne. 0) then
        call jeveuo(listOptionJv, 'L', vk16=listOption)
    end if

! - Get parameters for parallelism
    call getvtx(' ', 'PARALLELISME_TEMPS', scal=answer, nbret=iret)
    lTimeDist = answer .eq. 'OUI'

! - Prepare datastructure
    if (.not. lTimeDist) then
        call setResuPara(resultIn, resultOut, resultType, &
                         nbStore, listStore, listStoreJv, &
                         postComp%postCompResu)
    end if

! - Loop on options to compute
    do iOption = 1, nbOption
        option = listOption(iOption)
        call deprecated_option(option)
        if (option .eq. ' ') then
            cycle
        elseif ((option .eq. 'FORC_NODA') .or. (option .eq. 'REAC_NODA')) then
            if (lTimeDist) then
                call ccfnrn(option, resultIn, resultOut, listStoreJv, nbStore, &
                            resultType)
            else
                postComp%lPostNoda = ASTER_TRUE
                call ccfnrnLegacy(option, postComp)
            end if
        else
            call calcop(option, listOptionJv, resultIn, resultOut, listStoreJv, &
                        nbStore, resultType, iret, tldist_=.True._1)
            ASSERT(iret .eq. 0)
        end if
    end do
!
    call jedema()
!
end subroutine
