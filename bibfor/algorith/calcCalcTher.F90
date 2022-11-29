! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504
!
subroutine calcCalcTher(nb_option, list_option, &
                        list_load, model, mate, mateco, cara_elem, &
                        tpsthe, time, &
                        temp_prev, incr_temp, compor_ther, temp_curr, &
                        ve_charther, me_mtanther, ve_dirichlet, &
                        ve_evolther_l, ve_evolther_nl, ve_resither, &
                        nb_obje_maxi, obje_name, obje_sdname, nb_obje)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/knindi.h"
#include "asterfort/assert.h"
#include "asterfort/vechnl.h"
#include "asterfort/merxth.h"
#include "asterfort/vetnth_nonl.h"
#include "asterfort/verstp.h"
#include "asterfort/medith.h"
#include "asterfort/vethbt.h"
#include "asterfort/nmvcle.h"
!
    integer, intent(in) :: nb_option
    character(len=16), intent(in) :: list_option(:)
    character(len=19), intent(in) :: list_load
    character(len=24), intent(in) :: model, mate, mateco, cara_elem
    real(kind=8), intent(in) :: tpsthe(6)
    character(len=24), intent(in) :: time
    character(len=*), intent(in) :: temp_prev, incr_temp, temp_curr
    character(len=24), intent(in) :: compor_ther
    character(len=24), intent(in) :: ve_charther
    character(len=24), intent(inout) :: me_mtanther
    character(len=24), intent(in) :: ve_evolther_l, ve_evolther_nl
    character(len=*), intent(in) :: ve_dirichlet
    character(len=24), intent(inout) :: ve_resither
    integer, intent(in) :: nb_obje_maxi
    character(len=16), intent(inout) :: obje_name(nb_obje_maxi)
    character(len=24), intent(inout) :: obje_sdname(nb_obje_maxi)
    integer, intent(out) ::  nb_obje
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Compute for thermics
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_option        : number of options to compute
! In  list_option      : list of options to compute
! In  list_load        : name of datastructure for list of loads
! In  model            : name of model
! In  mate             : name of material characteristics (field)
! In  cara_elem        : name of elementary characteristics (field)
! In  tpsthe           : current time
! In  time             : name of field for time parameters
! In  temp_prev        : temperature at beginning of step
! In  incr_temp        : increment of temperature
! In  compor_ther      : name of comportment definition (field)
! In  temp_curr        : temperature at end of step
! In  ve_charther      : name of elementary for loads vector
! In  ve_mtanther      : name of elementary for tangent matrix
! In  ve_evolther_l    : name of elementary for linear transient vector
! In  ve_evolther_nl   : name of elementary for non-linear transient vector
! In  ve_resither      : name of elementary for residual vector
! In  ve_dirichlet     : name of elementary for reaction (Lagrange) vector
! In  nb_obje_maxi     : maximum number of new objects to add
! IO  obje_name        : name of new objects to add
! IO  obje_sdname      : datastructure name of new objects to add
! Out nb_obje          : number of new objects to add
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_char_ther, l_mtan_ther, l_char_evol, l_resi_ther, l_lagr, l_stat
    character(len=24) :: lload_name, lload_info, varc_curr
    real(kind=8) :: time_curr
!
! --------------------------------------------------------------------------------------------------
!
    lload_name = list_load(1:19)//'.LCHA'
    lload_info = list_load(1:19)//'.INFC'

    time_curr = tpsthe(1)
!
! - Prepare command variables
!
    varc_curr = '&&OP0026.VARC_CURR'
    call nmvcle(model, mate, cara_elem, time_curr, varc_curr)
!
! - What we are computing
!
    l_char_ther = (knindi(16, 'CHAR_THER_ELEM', list_option, nb_option) .gt. 0)
    l_mtan_ther = (knindi(16, 'MATR_TANG_ELEM', list_option, nb_option) .gt. 0)
    l_char_evol = (knindi(16, 'CHAR_EVOL_ELEM', list_option, nb_option) .gt. 0)
    l_resi_ther = (knindi(16, 'RESI_THER_ELEM', list_option, nb_option) .gt. 0)
    l_lagr = l_mtan_ther
!
! - Loads
!
    if (l_char_ther) then
        call vechnl(model, lload_name, lload_info, time, &
                    temp_prev, ve_charther, 'G')
    end if
!
! - Physical dof computation
!
    if (l_mtan_ther) then
        l_stat = ASTER_FALSE
        call merxth(model, lload_name, lload_info, cara_elem, mate, mateco, &
                    tpsthe, time, temp_curr, compor_ther, varc_curr, &
                    me_mtanther, 'G', l_stat)
    end if
!
! - Lagrange dof computation
!
    if (l_lagr) then
        call medith('G', 'CUMU', model, list_load, me_mtanther)
        call vethbt(model, lload_name, lload_info, cara_elem, mate, &
                    temp_curr, ve_dirichlet, 'G')
    end if
!
! - Loads
!
    if (l_char_evol) then
        call vetnth_nonl(model, cara_elem, mate, mateco, time, compor_ther, &
                         temp_prev, varc_curr, &
                         ve_evolther_l, ve_evolther_nl, 'G')
    end if
!
! - Residuals
!
    if (l_resi_ther) then
        call verstp(model, lload_name, lload_info, cara_elem, mateco, &
                    time_curr, time, compor_ther, temp_prev, incr_temp, &
                    varc_curr, ve_resither, 'G')
    end if
!
! - New objects in table
!
    nb_obje = 0
    if (l_lagr) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'FORC_DIRI_ELEM'
        obje_sdname(nb_obje) = ve_dirichlet
    end if
    if (l_char_ther) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'CHAR_THER_ELEM'
        obje_sdname(nb_obje) = ve_charther
    end if
    if (l_mtan_ther) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'MATR_TANG_ELEM'
        obje_sdname(nb_obje) = me_mtanther
    end if
    if (l_char_evol) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'CHAR_EVOL_ELEM'
        obje_sdname(nb_obje) = ve_evolther_l
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'CHAR_EVNL_ELEM'
        obje_sdname(nb_obje) = ve_evolther_nl
    end if
    if (l_resi_ther) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'RESI_THER_ELEM'
        obje_sdname(nb_obje) = ve_resither
    end if
!
end subroutine
