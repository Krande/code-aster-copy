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

subroutine liexco(sdcont, keywf, mesh, model, nb_cont_zone, &
                  nb_cont_elem, nb_cont_node)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/cfnbsf.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lireco.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: sdcont
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: nb_cont_zone
    integer(kind=8), intent(in) :: nb_cont_elem
    integer(kind=8), intent(in) :: nb_cont_node
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Save nodes and elements
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  keywf            : factor keyword to read
! In  mesh             : name of mesh
! In  model            : name of model
! In  nb_cont_zone     : number of zones of contact
! In  nb_cont_elem     : number of elements of contact
! In  nb_cont_node     : number of nodes of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdecma, jdecno
    integer(kind=8) :: nb_elem_slav, nb_node_slav
    integer(kind=8) :: nb_elem_mast, nb_node_mast
    integer(kind=8) :: i_zone, i_surf, i_elem, i_node
    character(len=24) :: list_elem_slav, list_elem_mast
    character(len=24) :: list_node_slav, list_node_mast
    integer(kind=8), pointer :: v_list(:) => null()
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=24) :: sdcont_noeuco
    integer(kind=8), pointer :: v_sdcont_noeuco(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    i_surf = 1
!
! - Temporary datastructures
!
    list_elem_mast = '&&LIEXCO.MAIL.MAIT'
    list_elem_slav = '&&LIEXCO.MAIL.ESCL'
    list_node_mast = '&&LIEXCO.NOEU.MAIT'
    list_node_slav = '&&LIEXCO.NOEU.ESCL'
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    sdcont_noeuco = sdcont_defi(1:16)//'.NOEUCO'
!
! - Create datastructure for nodes and elements
!
    call wkvect(sdcont_mailco, 'G V I', nb_cont_elem, vi=v_sdcont_mailco)
    call wkvect(sdcont_noeuco, 'G V I', nb_cont_node, vi=v_sdcont_noeuco)
!
! - Loop on zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Read
!
        call lireco(keywf, mesh, model, i_zone, list_elem_slav, &
                    list_elem_mast, list_node_slav, list_node_mast, nb_elem_slav, nb_node_slav, &
                    nb_elem_mast, nb_node_mast)
!
! ----- Master elements
!
        call cfnbsf(sdcont_defi, i_surf, 'MAIL', nb_elem_mast, jdecma)
        call jeveuo(list_elem_mast, 'L', vi=v_list)
        do i_elem = 1, nb_elem_mast
            v_sdcont_mailco(jdecma+i_elem) = v_list(i_elem)
        end do
!
! ----- Master nodes
!
        call cfnbsf(sdcont_defi, i_surf, 'NOEU', nb_node_mast, jdecno)
        call jeveuo(list_node_mast, 'L', vi=v_list)
        do i_node = 1, nb_node_mast
            v_sdcont_noeuco(jdecno+i_node) = v_list(i_node)
        end do
!
        i_surf = i_surf+1
!
! ----- Slave elements
!
        call cfnbsf(sdcont_defi, i_surf, 'MAIL', nb_elem_slav, jdecma)
        call jeveuo(list_elem_slav, 'L', vi=v_list)
        do i_elem = 1, nb_elem_slav
            v_sdcont_mailco(jdecma+i_elem) = v_list(i_elem)
        end do
!
! ----- Slave nodes
!
        call cfnbsf(sdcont_defi, i_surf, 'NOEU', nb_node_slav, jdecno)
        call jeveuo(list_node_slav, 'L', vi=v_list)
        do i_node = 1, nb_node_slav
            v_sdcont_noeuco(jdecno+i_node) = v_list(i_node)
        end do
!
        i_surf = i_surf+1
    end do
!
    call jedetr(list_elem_slav)
    call jedetr(list_elem_mast)
    call jedetr(list_node_slav)
    call jedetr(list_node_mast)

    call jedema()
!
end subroutine
