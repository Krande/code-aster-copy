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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgCreateCompIncr(compor, l_etat_init)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/copisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(inout) :: compor
    aster_logical, intent(in) :: l_etat_init
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Create CARTE Comportement for ELAS_INCR
!
! IO compor    : name of <CARTE> COMPORTEMENT
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_vale, nb_zone, nb_cmp_max, i_zone
    character(len=16) :: rela_comp, defo_comp, rela_new
    character(len=19) :: compor_elas
    integer(kind=8), pointer :: v_compor_desc(:) => null()
    character(len=16), pointer :: v_compor_vale(:) => null()
!
    call jemarq()
!
! --- Initialization
!
    compor_elas = "&&CGCOMP.COMPOR"
    call copisd('CHAMP_GD', 'V', compor, compor_elas)
!
! - Access to COMPOR
!
    call jeveuo(compor_elas//'.DESC', 'L', vi=v_compor_desc)
    call jeveuo(compor_elas//'.VALE', 'E', vk16=v_compor_vale)
    call jelira(compor_elas//'.VALE', 'LONMAX', nb_vale)
    nb_zone = v_compor_desc(3)
    nb_cmp_max = nb_vale/v_compor_desc(2)
!
    do i_zone = 1, nb_zone
        rela_comp = v_compor_vale(nb_cmp_max*(i_zone-1)+RELA_NAME)
        defo_comp = v_compor_vale(nb_cmp_max*(i_zone-1)+DEFO)
!
        if (rela_comp(1:4) .ne. "ELAS") then
            if (rela_comp(1:10) == "VMIS_ISOT_" .and. .not. (rela_comp(11:12) == "NL")) then
                rela_new = "ELAS_VMIS_"//rela_comp(11:16)
                v_compor_vale(nb_cmp_max*(i_zone-1)+RELA_NAME) = rela_new
                if (l_etat_init) then
                    ! v_compor_vale(nb_cmp_max*(i_zone-1)+INCRELAS)  = "COMP_INCR"
                    call utmess("F", "RUPTURE3_11")
                else
                    v_compor_vale(nb_cmp_max*(i_zone-1)+INCRELAS) = "COMP_ELAS"
                    v_compor_vale(nb_cmp_max*(i_zone-1)+NVAR) = "1"
                end if
            else
                call utmess("F", "RUPTURE3_8", sk=rela_comp)
            end if
        end if
        if ((defo_comp .ne. "PETIT") .and. (defo_comp .ne. "GREEN_LAGRANGE")) then
            call utmess("F", "RUPTURE3_9", sk=defo_comp)
        end if
    end do
!
    compor = compor_elas
!
    call jedema()
!
end subroutine
