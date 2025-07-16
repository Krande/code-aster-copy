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

subroutine dimecz(sdcont, mesh, nb_cont_zone)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnbsf.h"
#include "asterfort/cfzone.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mmnbnz.h"
#include "asterfort/mminfl.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: sdcont
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nb_cont_zone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Save contact counters - Counters by zone
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  mesh             : name of mesh
! In  nb_cont_zone     : number of zones of contact
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_methco
    integer(kind=8), pointer :: v_sdcont_methco(:) => null()
    integer(kind=8) :: zmeth
    integer(kind=8) :: i_zone
    integer(kind=8) :: jdecme, jdecmm, jdecne, jdecnm
    integer(kind=8) :: i_surf_escl, i_surf_mast
    integer(kind=8) :: nb_node_mast, nb_node_slav, nb_elem_mast, nb_elem_slav
    integer(kind=8) :: nb_node_mastc, nb_node_slavc, nb_elem_mastc, nb_elem_slavc
    integer(kind=8) :: nb_cont_poin
    character(len=24) :: sdcont_defi
    integer(kind=8) :: cont_form
    aster_logical :: l_verif
    integer(kind=8) :: nb_cont_poinc
!
! ----------------------------------------------------------------------
!
    nb_elem_slav = 0
    nb_node_slav = 0
    nb_elem_mast = 0
    nb_elem_slav = 0
    nb_cont_poin = 0
    nb_elem_slavc = 0
    nb_node_slavc = 0
    nb_elem_mastc = 0
    nb_elem_slavc = 0
    nb_cont_poinc = 0
    sdcont_defi = sdcont(1:8)//'.CONTACT'
!
! - Parameters
!
    cont_form = cfdisi(sdcont_defi, 'FORMULATION')
!
! - Datastructure for contact
!
    sdcont_methco = sdcont_defi(1:16)//'.METHCO'
    call jeveuo(sdcont_methco, 'E', vi=v_sdcont_methco)
    zmeth = cfmmvd('ZMETH')
!
! - Total number elements/zones and nodes / zones
!
    do i_zone = 1, nb_cont_zone
        call cfzone(sdcont_defi, i_zone, 'ESCL', i_surf_escl)
        call cfnbsf(sdcont_defi, i_surf_escl, 'MAIL', nb_elem_slav, jdecme)
        call cfnbsf(sdcont_defi, i_surf_escl, 'NOEU', nb_node_slav, jdecne)
        call cfzone(sdcont_defi, i_zone, 'MAIT', i_surf_mast)
        call cfnbsf(sdcont_defi, i_surf_mast, 'MAIL', nb_elem_mast, jdecmm)
        call cfnbsf(sdcont_defi, i_surf_mast, 'NOEU', nb_node_mast, jdecnm)
        l_verif = mminfl(sdcont_defi, 'VERIF', i_zone)
        if (l_verif) then
            nb_elem_slavc = 0
            nb_elem_mastc = 0
            nb_node_slavc = 0
            nb_node_mastc = 0
        else
            nb_elem_slavc = nb_elem_slav
            nb_elem_mastc = nb_elem_mast
            nb_node_slavc = nb_node_slav
            nb_node_mastc = nb_node_mast
        end if
        v_sdcont_methco(zmeth*(i_zone-1)+8) = nb_elem_slav
        v_sdcont_methco(zmeth*(i_zone-1)+9) = nb_node_slav
        v_sdcont_methco(zmeth*(i_zone-1)+10) = nb_elem_mast
        v_sdcont_methco(zmeth*(i_zone-1)+11) = nb_node_mast
        v_sdcont_methco(zmeth*(i_zone-1)+12) = nb_elem_slavc
        v_sdcont_methco(zmeth*(i_zone-1)+13) = nb_node_slavc
        v_sdcont_methco(zmeth*(i_zone-1)+14) = nb_elem_mastc
        v_sdcont_methco(zmeth*(i_zone-1)+15) = nb_node_mastc
    end do
!
! - Shift/zone
!
    do i_zone = 1, nb_cont_zone
        call cfzone(sdcont_defi, i_zone, 'ESCL', i_surf_escl)
        call cfnbsf(sdcont_defi, i_surf_escl, 'MAIL', nb_elem_slav, jdecme)
        call cfnbsf(sdcont_defi, i_surf_escl, 'NOEU', nb_node_slav, jdecne)
        call cfzone(sdcont_defi, i_zone, 'MAIT', i_surf_mast)
        call cfnbsf(sdcont_defi, i_surf_mast, 'MAIL', nb_elem_mast, jdecmm)
        call cfnbsf(sdcont_defi, i_surf_mast, 'NOEU', nb_node_mast, jdecnm)
        v_sdcont_methco(zmeth*(i_zone-1)+16) = jdecme
        v_sdcont_methco(zmeth*(i_zone-1)+17) = jdecmm
        v_sdcont_methco(zmeth*(i_zone-1)+18) = jdecne
        v_sdcont_methco(zmeth*(i_zone-1)+19) = jdecnm
    end do
!
! - Number of contact points/zone
!
    do i_zone = 1, nb_cont_zone
        call mmnbnz(mesh, sdcont_defi, i_zone, nb_cont_poin)
        v_sdcont_methco(zmeth*(i_zone-1)+20) = nb_cont_poin
    end do
!
! - No computation mode
!
    do i_zone = 1, nb_cont_zone
        nb_cont_poin = v_sdcont_methco(zmeth*(i_zone-1)+20)
        l_verif = mminfl(sdcont_defi, 'VERIF', i_zone)
        if (l_verif) then
            nb_cont_poinc = 0
        else
            nb_cont_poinc = nb_cont_poin
        end if
        v_sdcont_methco(zmeth*(i_zone-1)+21) = nb_cont_poinc
    end do
!
end subroutine
