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
subroutine nonlinDSInOutRead(phenom, result, ds_inout)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/gcucon.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    character(len=4), intent(in) :: phenom
    character(len=8), intent(in) :: result
    type(NL_DS_InOut), intent(inout) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE and THER_* - Input/output management
!
! Read parameters for input/output management
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom           : phenomenon (MECA/THER/ACOU)
! In  result           : name of result datastructure (EVOL_NOLI)
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=16) :: keywf, init_keyw, result_type
    character(len=8) :: stin_evol, field, criterion, answer
    real(kind=8) :: precision, user_time, temp_init
    integer(kind=8) :: nocc, didi_nume, i_field, user_nume, iret
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_3')
    end if
!
! - Initializations
!
    keywf = 'ETAT_INIT'
    if (phenom .eq. 'MECA') then
        result_type = 'EVOL_NOLI'
    elseif (phenom .eq. 'THER') then
        result_type = 'EVOL_THER'
    elseif (phenom .eq. 'VIBR') then
        result_type = 'DYNA_TRANS'
    else
        ASSERT(.false.)
    end if
!
! - Set name of result datastructure
!
    ds_inout%result = result
    if (getexm(keywf, result_type) .ne. 1) then
        goto 99
    end if
!
! - Does ETAT_INIT (initial state) exist ?
!
    call getfac(keywf, nocc)
    ASSERT(nocc .le. 1)
    ds_inout%l_state_init = nocc .gt. 0
!
! - Is REUSE ?
!
    call gcucon(result, result_type, iret)
    ds_inout%l_reuse = iret .gt. 0
!
! - Get name of result datastructure in ETAT_INIT
!
    call getvid(keywf, result_type, iocc=1, scal=stin_evol, nbret=nocc)
    if (nocc .ge. 1) then
        ds_inout%stin_evol = stin_evol
        ds_inout%l_stin_evol = .true.
    end if

    if (ds_inout%l_reuse .and. (.not. ds_inout%l_state_init)) then
        call utmess('F', 'MECANONLINE_4')
    end if
!
! - For DIDI loads
!
    if (phenom .eq. 'MECA') then
        call getvis(keywf, 'NUME_DIDI', iocc=1, scal=didi_nume, nbret=nocc)
        if (nocc .ge. 1) then
            ds_inout%didi_nume = didi_nume
        end if
    end if
!
! - For thermics
!
    if (phenom .eq. 'THER') then
        call getvtx(keywf, 'STAT', iocc=1, scal=answer, nbret=nocc)
        if (nocc .eq. 1) then
            ds_inout%l_init_stat = ASTER_TRUE
        end if
        call getvr8(keywf, 'VALE', iocc=1, scal=temp_init, nbret=nocc)
        if (nocc .eq. 1) then
            ds_inout%l_init_vale = ASTER_TRUE
            ds_inout%temp_init = temp_init
        end if
    end if
!
! - Read fields
!
    do i_field = 1, ds_inout%nb_field
        init_keyw = ds_inout%field(i_field)%init_keyw
        if (getexm(keywf, init_keyw) .eq. 1 .and. ds_inout%field(i_field)%l_read_init) then
            call getvid(keywf, init_keyw, scal=field, iocc=1, nbret=nocc)
            if (nocc .eq. 1) then
                ds_inout%l_field_read(i_field) = ASTER_TRUE
                ds_inout%field(i_field)%field_read = field
            end if
        end if
    end do
!
! - Get parameters to read initial state
!
    call getvr8(keywf, 'PRECISION', iocc=1, scal=precision, nbret=nocc)
    if (nocc .eq. 0) then
        ds_inout%precision = 1.d-6
    else
        ds_inout%precision = precision
    end if
    call getvtx(keywf, 'CRITERE', iocc=1, scal=criterion, nbret=nocc)
    if (nocc .eq. 0) then
        ds_inout%criterion = 'RELATIF'
    else
        ds_inout%criterion = criterion
    end if

    call getvis(keywf, 'NUME_ORDRE', iocc=1, scal=user_nume, nbret=nocc)
    ds_inout%user_nume = user_nume
    ds_inout%l_user_nume = nocc .gt. 0

    if (phenom .eq. 'VIBR') then
        call getvr8(keywf, 'INST_INIT', iocc=1, scal=user_time, nbret=nocc)
        ds_inout%user_time = user_time
        ds_inout%l_user_time = nocc .gt. 0
        ds_inout%stin_time = user_time
        ds_inout%l_stin_time = nocc .gt. 0
    else
        ! MECA or THER
        call getvr8(keywf, 'INST', iocc=1, scal=user_time, nbret=nocc)
        ds_inout%user_time = user_time
        ds_inout%l_user_time = nocc .gt. 0
        ds_inout%stin_time = 0.d0
        ds_inout%l_stin_time = ASTER_FALSE
    end if
!
99  continue
!
end subroutine
