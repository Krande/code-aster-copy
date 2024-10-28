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

subroutine load_neum_evcm(inst_curr, load_name, i_load, ligrel_calc, &
                          nb_in_maxi, nb_in_prep, lpain, lchin, &
                          idx_matr, matr_elem)
!
    implicit none
!
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/calcul.h"
#include "asterfort/chpnua.h"
#include "asterfort/cnocre.h"
#include "asterfort/copisd.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcncon.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/nuachp.h"
#include "asterfort/pronua.h"
#include "asterfort/rsinch.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtgpld.h"
#include "asterfort/load_neum_matr.h"
!
!
    character(len=24), intent(in)  :: ligrel_calc
    real(kind=8), intent(in) :: inst_curr
    character(len=8), intent(in) :: load_name
    integer, intent(in) :: i_load
    integer, intent(in) :: nb_in_maxi
    character(len=*), intent(inout) :: lpain(nb_in_maxi)
    character(len=*), intent(inout) :: lchin(nb_in_maxi)
    integer, intent(in) :: nb_in_prep
    integer, intent(inout) :: idx_matr
    character(len=19), intent(in) :: matr_elem
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann matrix
!
! EVOL_CHAR - Undead loads
!
! --------------------------------------------------------------------------------------------------
!
! In  ligrel_calc    : LIGREL to compute
! In  inst_curr      : current time
! In  i_load         : index of current load
! In  nb_in_maxi     : maximum number of input fields
! In  nb_in_prep     : number of input fields before specific ones
! IO  lpain          : list of input parameters
! IO  lchin          : list of input fields
! In  load_name      : name of current load
! In  matr_elem      : name of vect_elem
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ier, nb_cham
    character(len=4) :: load_type
    character(len=8) :: evol_char
    character(len=16) :: type_sd, option
    character(len=19) :: load_name_evol, iden_direct
    integer :: load_nume_evol
    character(len=24) :: object
    character(len=8), pointer :: p_object(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    load_name_evol = '&&NMMGME'
!
! - Only scalar loadings, dead and fixed loads
!
    load_nume_evol = 1
    load_type = 'Suiv'
!
! - Get evol_char
!
    object = load_name//'.CHME.EVOL.CHAR'
    call jeexin(object, ier)
    if (ier .eq. 0) then
        goto 99
    end if
    call jeveuo(object, 'L', vk8=p_object)
    evol_char = p_object(1)
!
! - Check
!
    call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nb_cham)
    ASSERT(nb_cham .gt. 0)
    call gettco(evol_char, type_sd)
    ASSERT(type_sd .eq. 'EVOL_CHAR')
!
! - Get pressure (CHAR_MECA_PRES_R)
!
    call rsinch(evol_char, 'PRES', 'INST', inst_curr, load_name_evol, &
                'EXCLU', 'EXCLU', 0, 'V', ier)
    if (ier .le. 2) then
        option = 'RIGI_MECA_PRSU_R'
        goto 30
    else if (ier .eq. 11 .or. ier .eq. 12 .or. ier .eq. 20) then
        call utmess('F', 'CHARGES3_8', sk=evol_char, sr=inst_curr)
    end if
30  continue
!
! - Compute pressure (RIGI_MECA_PRSU_R)
!
    if (option .eq. 'RIGI_MECA_PRSU_R') then
        iden_direct = '.PRESS'
        call load_neum_matr(i_load, idx_matr, load_name, load_nume_evol, load_type, &
                            ligrel_calc, nb_in_maxi, nb_in_prep, lpain, lchin, &
                            matr_elem, iden_direct=iden_direct, &
                            name_inputz=load_name_evol)
    end if
!
99  continue
!
    call jedema()
end subroutine
