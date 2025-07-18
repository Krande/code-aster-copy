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
! aslint: disable=W0413
!
subroutine cazofm(sdcont, keywf, cont_form, cont_nbzone)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cazouu.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: sdcont
    integer(kind=8), intent(in) :: cont_form
    integer(kind=8), intent(in) :: cont_nbzone
    character(len=16), intent(in) :: keywf
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get method of contact
!
! --------------------------------------------------------------------------------------------------
!
! In  keywf       : factor keyword to read
! In  sdcont      : name of contact concept (DEFI_CONTACT)
! In  cont_form   : formulation of contact
! In  cont_nbzone : number of zones of contact
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: s_algo_cont, s_formul, s_frott
    integer(kind=8) :: algo_cont, algo_frot, izone
    aster_logical :: l_frot, lmunul
    character(len=24) :: sdcont_paraci
    integer(kind=8), pointer :: v_sdcont_paraci(:) => null()
    real(kind=8) :: coefff
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    algo_cont = 0
    algo_frot = 0
    l_frot = ASTER_FALSE
    s_formul = ' '
    s_frott = ' '
    s_algo_cont = ' '
    if (cont_form .eq. 1) then
        s_formul = 'DISCRETE'
    else if (cont_form .eq. 2) then
        s_formul = 'CONTINUE'
    else if (cont_form .eq. 5) then
        s_formul = 'LAC'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Datastructure for contact definition
    sdcont_paraci = sdcont(1:8)//'.PARACI'
    call jeveuo(sdcont_paraci, 'E', vi=v_sdcont_paraci)
!
! - Formulation
!
    v_sdcont_paraci(4) = cont_form

! - Friction ?
    call getvtx(' ', 'FROTTEMENT', scal=s_frott)
    l_frot = s_frott .eq. 'COULOMB'

! - Get methods
    call getvtx(keywf, 'ALGO_CONT', iocc=1, scal=s_algo_cont)
    if (cont_form .eq. 1) then
        call cazouu(keywf, cont_nbzone, 'ALGO_CONT', 'T')
        if (l_frot) then
            algo_frot = 1
            if (s_algo_cont .eq. 'PENALISATION') then
                algo_cont = 4
            else
                call utmess('F', 'CONTACT3_3')
            end if
        else
            algo_frot = 0
            if (s_algo_cont .eq. 'GCP') then
                algo_cont = 2
            else if (s_algo_cont .eq. 'CONTRAINTE') then
                algo_cont = 1
            else if (s_algo_cont .eq. 'PENALISATION') then
                algo_cont = 4
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    else if (cont_form .eq. 2) then
        algo_cont = 6
        if (l_frot) then
            lmunul = ASTER_FALSE
            do izone = 1, cont_nbzone
                call getvr8(keywf, 'COULOMB', iocc=izone, scal=coefff)
                lmunul = lmunul .or. (coefff .ne. 0.d0)
            end do
            if (.not. lmunul) then
                call utmess('A', 'CONTACT3_1')
                l_frot = ASTER_FALSE
            end if
        end if
        if (l_frot) then
            algo_frot = 6
        else
            algo_frot = 0
        end if
    else if (cont_form .eq. 5) then
        algo_cont = 8
    else
        ASSERT(ASTER_FALSE)
    end if

! - Save methods
    v_sdcont_paraci(17) = algo_cont
    v_sdcont_paraci(18) = algo_frot
!
    call jedema()
!
end subroutine
