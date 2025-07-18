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

subroutine ntcrob(meshz, modelz, result, sddisc, ds_inout, &
                  sd_obsv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcroi.h"
#include "asterfort/nmcrot.h"
#include "asterfort/nmextr.h"
#include "asterfort/nmobno.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=*), intent(in) :: meshz
    character(len=*), intent(in) :: modelz
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: sddisc
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(out) :: sd_obsv
!
! --------------------------------------------------------------------------------------------------
!
! THER_* - Observation
!
! Create observation datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  result           : name of results datastructure
! In  sddisc           : datastructure for discretization
! In  ds_inout         : datastructure for input/output management
! Out sd_obsv          : datastructure for observation parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_obsv, nb_keyw_fact, numeReuseCalc
    character(len=19) :: sdarch
    character(len=14) :: sdextr_obsv
    character(len=16) :: keyw_fact
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    character(len=24) :: extr_info
    integer(kind=8), pointer :: v_extr_info(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_obsv = 0
    sd_obsv = '&&NMCROB.OBSV'
    keyw_fact = 'OBSERVATION'
    call getfac(keyw_fact, nb_keyw_fact)
    ASSERT(nb_keyw_fact .le. 99)
!
! - Access to storage datastructure
!
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    call jeveuo(sdarchAinfJv, 'L', vi=sdarchAinf)
!
! - Read datas for extraction
!
    sdextr_obsv = sd_obsv(1:14)
    call nmextr(meshz, modelz, sdextr_obsv, ds_inout, keyw_fact, &
                nb_keyw_fact, nb_obsv)
!
! - Set reuse index in OBSERVATION table
!
    numeReuseCalc = sdarchAinf(3)
    extr_info = sdextr_obsv(1:14)//'     .INFO'
    call jeveuo(extr_info, 'E', vi=v_extr_info)
    v_extr_info(4) = numeReuseCalc
!
! - Read parameters
!
    if (nb_obsv .ne. 0) then
!
        call utmess('I', 'OBSERVATION_3', si=nb_obsv)
!
! ----- Read time list
!
        call nmcroi(sd_obsv, keyw_fact, nb_keyw_fact)
!
! ----- Read name of columns
!
        call nmobno(sd_obsv, keyw_fact, nb_keyw_fact)
!
! ----- Create table
!
        call nmcrot(result, sd_obsv)
    end if
!
end subroutine
