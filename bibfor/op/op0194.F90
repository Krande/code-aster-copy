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
subroutine op0194()
!
    use MetallurgyOperator_module
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/detrsd.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mtdorc.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rslesd.h"
#include "asterfort/smevol.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! Command: CALC_META
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iret, n1, numeFieldInit
    integer :: nbOption, nb, iOption
    integer :: nbStore, numeStore0
    character(len=8) :: resultName, model, materialField
    character(len=16) :: resultType, option
    character(len=24), parameter :: comporMeta = '&&OP0194.COMPOR'
    character(len=24) :: metaInitUser, materialCoding
    character(len=24), parameter :: listOptionsJv = '&&OP0194.LES_OPTION'
    character(len=16), parameter :: keywordfact = 'COMPORTEMENT'
    integer :: i_comp, nbocc
    character(len=16) :: phase_type
    character(len=19), parameter :: listStoreJv = '&&OP0194.LISTSTORE'
    integer, pointer :: listStore(:) => null()
    character(len=16), pointer :: listOption(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Get output result
    call getvid(' ', 'RESULTAT', scal=resultName, nbret=n1)
    call gettco(resultName, resultType)
    ASSERT(resultType .eq. 'EVOL_THER')

! - Get all storing index in result
    call rs_get_liststore(resultName, nbStore)
    if (nbStore .ne. 0) then
        call wkvect(listStoreJv, "V V I", nbStore, vi=listStore)
        call rs_get_liststore(resultName, nbStore, listStore)
    end if
    if (nbStore .lt. 2) then
        call utmess('F', 'META1_1')
    end if

! - Get main parameters
    materialCoding = ' '
    model = ' '
    materialField = ' '
    numeStore0 = listStore(1)
    call rslesd(resultName, numeStore0, model, materialField)
    if (materialField .ne. ' ') then
        call rcmfmc(materialField, materialCoding, l_ther_=ASTER_TRUE)
    end if
!
! - Get options to compute
    call getvtx(' ', 'OPTION', nbval=0, nbret=nb)
    nbOption = -nb
    call wkvect(listOptionsJv, "V V K16", nbOption, vk16=listOption)
    call getvtx(' ', 'OPTION', nbval=nbOption, vect=listOption, nbret=nb)

! - Compute options
    do iOption = 1, nbOption
        option = listOption(iOption)
!
        if (option .eq. 'META_ELNO') then

! --------- Construct map for thermic behaviour
            call mtdorc(model, comporMeta)

! --------- Initial state
            call metaGetInitialState(resultName, metaInitUser, numeFieldInit)

! --------- Compute
            call smevol(resultName, nbStore, listStore, &
                        model, materialField, materialCoding, comporMeta, &
                        metaInitUser, numeFieldInit)
!
            call detrsd('CARTE', comporMeta)
!
        else
            nbocc = 0
            call getfac(keywordfact, nbocc)
            do i_comp = 1, nbocc
                call getvtx(keywordfact, 'RELATION', iocc=i_comp, scal=phase_type, nbret=iret)
                if (iret .ne. 0) then
                    if (phase_type(1:5) .ne. 'ACIER' .and. option .eq. 'DURT_ELNO') then
                        call utmess('F', 'META1_3', sk=phase_type)
                    end if
                end if
            end do
            call calcop(option, listOptionsJv, resultName, resultName, listStoreJv, &
                        nbStore, resultType, iret)
            if (iret .eq. 0) cycle
        end if
    end do

! - Cleaning
    call jedetr(listOptionsJv)
    call jedetr(listStoreJv)
!
    call jedema()
end subroutine
