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

subroutine surfll(sdcont_defi, mesh, unit_msg, nb_cont_zone, nb_cont_elem, &
                  nb_cont_node)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfnbsf.h"
#include "asterfort/cfzone.h"
#include "asterfort/jeveuo.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: sdcont_defi
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: unit_msg
    integer(kind=8), intent(in) :: nb_cont_elem
    integer(kind=8), intent(in) :: nb_cont_node
    integer(kind=8), intent(in) :: nb_cont_zone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Print list of elements and list of nodes
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  mesh             : name of mesh
! In  unit_msg         : logical unit for messages (print)
! In  nb_cont_zone     : number of zones of contact
! In  nb_cont_node     : number of nodes of contact
! In  nb_cont_elem     : number of elements of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_elem, nb_node
    integer(kind=8) :: i_zone, i_surf, i_elem, i_node
    integer(kind=8) :: jdecma, jdecno
    integer(kind=8) :: elem_nume, node_nume
    character(len=8) :: chain1, chain2
    character(len=24) :: sdcont_pzoneco
    integer(kind=8), pointer :: v_sdcont_pzoneco(:) => null()
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=24) :: sdcont_noeuco
    integer(kind=8), pointer :: v_sdcont_noeuco(:) => null()
    character(len=24) :: sdcont_psumaco
    integer(kind=8), pointer :: v_sdcont_psumaco(:) => null()
    character(len=24) :: sdcont_psunoco
    integer(kind=8), pointer :: v_sdcont_psunoco(:) => null()
    character(len=8), pointer :: v_trav_elem(:) => null()
    character(len=8), pointer :: v_work_node(:) => null()
!
! ----------------------------------------------------------------------
!
!
! - Datastructure for contact definition
!
    sdcont_pzoneco = sdcont_defi(1:16)//'.PZONECO'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    sdcont_noeuco = sdcont_defi(1:16)//'.NOEUCO'
    sdcont_psumaco = sdcont_defi(1:16)//'.PSUMACO'
    sdcont_psunoco = sdcont_defi(1:16)//'.PSUNOCO'
    call jeveuo(sdcont_pzoneco, 'L', vi=v_sdcont_pzoneco)
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
    call jeveuo(sdcont_noeuco, 'L', vi=v_sdcont_noeuco)
    call jeveuo(sdcont_psumaco, 'L', vi=v_sdcont_psumaco)
    call jeveuo(sdcont_psunoco, 'L', vi=v_sdcont_psunoco)
!
! - Temporary vectors
!
    AS_ALLOCATE(vk8=v_trav_elem, size=nb_cont_elem)
    AS_ALLOCATE(vk8=v_work_node, size=nb_cont_node)
!
! - Global parameters
!
    write (unit_msg, *)
    write (unit_msg, *) '<CONTACT> INFOS SUR LES SURFACES MAILLEES '
    write (unit_msg, *)
    do i_zone = 1, nb_cont_zone
!
! ----- Zone
!
        write (unit_msg, *) '<CONTACT> ZONE : ', i_zone
        nb_elem = v_sdcont_psumaco(v_sdcont_pzoneco(i_zone+1)+1)- &
                  v_sdcont_psumaco(v_sdcont_pzoneco(i_zone)+1)
        write (unit_msg, *) '<CONTACT> ... NOMBRE DE MAILLES          : ', nb_elem
        nb_node = v_sdcont_psunoco(v_sdcont_pzoneco(i_zone+1)+1)- &
                  v_sdcont_psunoco(v_sdcont_pzoneco(i_zone)+1)
        write (unit_msg, *) '<CONTACT> ... NOMBRE DE NOEUDS           : ', nb_node
!
! ----- Master surface
!
        write (unit_msg, *) '<CONTACT> ...... SURFACE MAITRE '
        call cfzone(sdcont_defi, i_zone, 'MAIT', i_surf)
        call cfnbsf(sdcont_defi, i_surf, 'MAIL', nb_elem, jdecma)
        call cfnbsf(sdcont_defi, i_surf, 'NOEU', nb_node, jdecno)
        if (nb_elem .le. 1) then
            chain1 = ' MAILLE '
        else
            chain1 = ' MAILLES'
        end if
        if (nb_node .le. 1) then
            chain2 = ' NOEUD  '
        else
            chain2 = ' NOEUDS '
        end if
        write (unit_msg, 135) nb_elem, chain1, ' ET ', nb_node, chain2
        do i_elem = 1, nb_elem
            elem_nume = v_sdcont_mailco(jdecma+i_elem)
            v_trav_elem(i_elem) = int_to_char8(elem_nume)
        end do
        write (unit_msg, 104) '     LISTE DES MAILLES : '
        write (unit_msg, 105) (v_trav_elem(i_elem), i_elem=1, nb_elem)
        do i_node = 1, nb_node
            node_nume = v_sdcont_noeuco(jdecno+i_node)
            v_work_node(i_node) = int_to_char8(node_nume)
        end do
        write (unit_msg, 104) '     LISTE DES NOEUDS  : '
        write (unit_msg, 105) (v_work_node(i_node), i_node=1, nb_node)
!
! ----- Slave surface
!
        write (unit_msg, *) '<CONTACT> ...... SURFACE ESCLAVE '
        call cfzone(sdcont_defi, i_zone, 'ESCL', i_surf)
        call cfnbsf(sdcont_defi, i_surf, 'MAIL', nb_elem, jdecma)
        call cfnbsf(sdcont_defi, i_surf, 'NOEU', nb_node, jdecno)
        if (nb_elem .le. 1) then
            chain1 = ' MAILLE '
        else
            chain1 = ' MAILLES'
        end if
        if (nb_node .le. 1) then
            chain2 = ' NOEUD  '
        else
            chain2 = ' NOEUDS '
        end if
        write (unit_msg, 135) nb_elem, chain1, ' ET ', nb_node, chain2
        do i_elem = 1, nb_elem
            elem_nume = v_sdcont_mailco(jdecma+i_elem)
            v_trav_elem(i_elem) = int_to_char8(elem_nume)
        end do
        write (unit_msg, 104) '     LISTE DES MAILLES : '
        write (unit_msg, 105) (v_trav_elem(i_elem), i_elem=1, nb_elem)
        do i_node = 1, nb_node
            node_nume = v_sdcont_noeuco(jdecno+i_node)
            v_work_node(i_node) = int_to_char8(node_nume)
        end do
        write (unit_msg, 104) '     LISTE DES NOEUDS  : '
        write (unit_msg, 105) (v_work_node(i_node), i_node=1, nb_node)
    end do
!
135 format(' <CONTACT> ...... ', i5, a8, a4, i5, a8)
104 format(' <CONTACT> ...... ', a25)
105 format((' <CONTACT> ...... ', 17x, 4(a8, 1x)))
!
! - Clean
!
    AS_DEALLOCATE(vk8=v_trav_elem)
    AS_DEALLOCATE(vk8=v_work_node)
!
end subroutine
