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
subroutine op0026()
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/catabl.h"
#include "asterfort/infmaj.h"
#include "asterfort/jemarq.h"
#include "asterfort/diinst.h"
#include "asterfort/nmchai.h"
#include "asterfort/jedema.h"
#include "asterfort/calcGetData.h"
#include "asterfort/calcGetDataMeca.h"
#include "asterfort/calcPrepDataMeca.h"
#include "asterfort/calcCalcMeca.h"
!
! --------------------------------------------------------------------------------------------------
!
!  O P E R A T E U R    C A L C U L
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: zsolal = 17, zvalin = 28
    character(len=19) :: hval_incr(zvalin), hval_algo(zsolal)
    integer(kind=8) :: nb_obje
    integer(kind=8), parameter :: nb_obje_maxi = 9
    character(len=16) :: obje_name(nb_obje_maxi)
    character(len=24) :: obje_sdname(nb_obje_maxi)
    integer(kind=8) :: nb_option
    character(len=16) :: list_option(6)
    integer(kind=8) :: nume_inst, nume_harm
    integer(kind=8) :: long
    real(kind=8) :: time_prev, time_curr
    real(kind=8) :: deltat, khi
    character(len=8) :: table_new, table_old
    type(NL_DS_Constitutive) :: ds_constitutive
    character(len=16) :: phenom
    character(len=19) :: list_inst
    character(len=8) :: model, materField, caraElem
    character(len=24) :: mateco, listLoad
    character(len=19) :: disp_prev, disp_cumu_inst, vari_prev, sigm_prev
    character(len=19) :: vediri, vefnod, vevarc_prev, vevarc_curr
    aster_logical :: l_elem_nonl
    type(NL_DS_Material) :: ds_material
    type(NL_DS_System) :: ds_system
    aster_logical :: l_pred

! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
! - Get commons data
!
    call calcGetData(table_new, table_old, &
                     nb_option, list_option, &
                     nume_inst, list_inst, &
                     phenom, l_pred)
    if (phenom .eq. 'MECANIQUE') then
        call calcGetDataMeca(listLoad, model, materField, mateco, caraElem, &
                             disp_prev, disp_cumu_inst, vari_prev, sigm_prev, &
                             ds_constitutive, l_elem_nonl, nume_harm)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Get current and previous times
!
    time_prev = diinst(list_inst, nume_inst-1)
    time_curr = diinst(list_inst, nume_inst)
    deltat = time_curr-time_prev
    khi = 1.d0
!
! - Check lengths
!
    call nmchai('VALINC', 'LONMAX', long)
    ASSERT(long .eq. zvalin)
    call nmchai('SOLALG', 'LONMAX', long)
    ASSERT(long .eq. zsolal)
!
! - Prepare data
!
    if (phenom .eq. 'MECANIQUE') then
        call calcPrepDataMeca(model, materField, mateco, caraElem, &
                              disp_prev, disp_cumu_inst, vari_prev, sigm_prev, &
                              time_prev, time_curr, &
                              ds_constitutive, ds_material, ds_system, &
                              hval_incr, hval_algo, &
                              vediri, vefnod, &
                              vevarc_prev, vevarc_curr)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Compute
!
    if (phenom .eq. 'MECANIQUE') then
        call calcCalcMeca(nb_option, list_option, &
                          l_elem_nonl, &
                          listLoad, model, caraElem, &
                          ds_constitutive, ds_material, ds_system, &
                          hval_incr, hval_algo, &
                          vediri, vefnod, &
                          nb_obje_maxi, obje_name, obje_sdname, nb_obje, &
                          l_pred)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Table management
!
    if (phenom .eq. 'MECANIQUE') then
        call catabl(table_new, table_old, time_curr, nume_inst, nb_obje, &
                    obje_name, obje_sdname)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
!
end subroutine
