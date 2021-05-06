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
!
subroutine op0194()
!
implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/modopt.h"
#include "asterfort/mtdorc.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rsexch.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rs_getnume.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsorac.h"
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
    integer :: iret, n1, n2, n3, numpha
    integer :: nbtrou, ier,  nbOption, nb, iOption
    integer :: nbStore, numeStore0, numeStoreInit
    real(kind=8) :: inst, prec
    character(len=4) :: ctyp
    character(len=8) :: crit, temper, temper2, model
    character(len=16) :: resultType, option
    character(len=19), parameter :: compor = '&&OP0194.COMPOR'
    character(len=24) :: chmeta, phasin, mateco, chmate
    character(len=24), parameter :: listOptionsJv = '&&OP0194.LES_OPTION'
    character(len=16), parameter :: keywordfact  = 'COMPORTEMENT'
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
    call getvid(' ', 'RESULTAT', scal=temper, nbret=n1)
    call gettco(temper, resultType)
    ASSERT(resultType .eq. 'EVOL_THER')
    ctyp = ' '

! - Get list of storing index
    call rs_get_liststore(temper, nbStore)
    if (nbStore .lt. 2) then
        call utmess('F', 'META1_1')
    endif
    call wkvect(listStoreJv, 'V V I', nbStore, vi = listStore)
    call rs_get_liststore(temper, nbStore, listStore)
    numeStore0 = listStore(1)

! - Get main parameters
    mateco = ' '
    model  = ' '
    chmate = ' '
    call rslesd(temper, numeStore0, model, chmate)
    if (chmate .ne. ' ') then
        call rcmfmc(chmate, mateco, l_ther_ = ASTER_TRUE)
    endif
!
! - Get options to compute
    call getvtx(' ', 'OPTION', nbval=0, nbret=nb)
    nbOption = -nb
    AS_ALLOCATE(vk16 = listOption, size = nbOption)
    call getvtx(' ', 'OPTION', nbval=nbOption, vect=listOption, nbret=nb)

! - Compute options
    do iOption = 1, nbOption
        option = listOption(iOption)
!
        if (option .eq. 'META_ELNO') then

! --------- Construct map for thermic behaviour
            call mtdorc(model, compor)

! --------- Initial state
            numpha = 0
            call getvid('ETAT_INIT', 'META_INIT_ELNO', iocc=1, scal=chmeta, nbret=n3)
            if (n3 .gt. 0) then
                phasin = '&&SMEVOL_ZINIT'
                call chpver('F', chmeta(1:19), 'CART', 'VAR2_R', ier)
                call copisd('CHAMP_GD', 'V', chmeta, phasin)
            else
                call getvid('ETAT_INIT', 'EVOL_THER', iocc=1, scal=temper2, nbret=n1)
                if (temper2 .ne. temper) then
                    call utmess('F', 'META1_2')
                endif
                call getvis('ETAT_INIT', 'NUME_INIT', iocc=1, scal=numeStoreInit, nbret=n2)
                if (n2 .eq. 0) then
                    call getvr8('ETAT_INIT', 'INST_INIT', iocc=1, scal=inst, nbret=n3)
                    call getvr8('ETAT_INIT', 'PRECISION', iocc=1, scal=prec, nbret=n3)
                    call getvtx('ETAT_INIT', 'CRITERE', iocc=1, scal=crit, nbret=n3)
                    call rs_getnume(temper, inst, crit, prec, numeStoreInit, nbtrou)
                    if (nbtrou .eq. 0) then
                        call utmess('F', 'UTILITAI6_51', sk=temper, sr=inst)
                    else if (nbtrou .gt. 1) then
                        call utmess('F', 'UTILITAI6_52', sk=temper, si=nbtrou, sr=inst)
                    endif
                endif
                call rsexch('F', temper, 'META_ELNO', numeStoreInit, phasin, iret)
                numpha = numeStoreInit
            endif

! --------- Compute
            call smevol(temper, model, chmate, mateco, compor, option,&
                        phasin, numpha)
!
            call detrsd('CARTE', '&&NMDORC.COMPOR')
!
        else
            nbocc = 0
            call getfac(keywordfact, nbocc)
            do i_comp = 1, nbocc
                call getvtx(keywordfact, 'RELATION', iocc = i_comp, scal = phase_type, nbret=iret)
                if (iret .ne. 0) then
                    if (phase_type(1:5) .ne. 'ACIER' .and. option.eq.'DURT_ELNO') then
                        call utmess('F', 'META1_3', sk=phase_type)
                    endif
                endif
            end do
            call calcop(option, listOptionsJv, temper, temper, listStoreJv,&
                        nbStore, ctyp, resultType, iret)
            if (iret .eq. 0) cycle
!
        endif
    end do
!
    AS_DEALLOCATE(vk16 = listOption)
!
    call jedema()
end subroutine
