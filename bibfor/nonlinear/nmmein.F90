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

subroutine nmmein(mesh, model, crack, nb_dim, list_node, &
                  nb_node, list_cmp, list_node_1, list_node_2, cmp_name, &
                  nb_node_sele)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmaret.h"
#include "asterfort/wkvect.h"
#include "asterfort/xlagsp.h"
!
!
    integer, intent(in) :: nb_dim
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    character(len=8), intent(in)  :: crack
    integer, intent(in) :: nb_node
    character(len=24), intent(in) :: list_node
    character(len=24), intent(in) :: list_cmp
    character(len=24), intent(in) :: list_node_1
    character(len=24), intent(in) :: list_node_2
    character(len=8), intent(out) :: cmp_name
    integer, intent(out) :: nb_node_sele
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Initializations
!
! Select edges and component for continuation method in XFEM
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh           : name of mesh
! In  model          : name of model
! In  crack          : name of crack
! In  nb_dim         : dimension of space
! In  nb_node        : number of nodes in list_node
! In  list_node      : list of nodes to apply continuation
! In  list_cmp       : list of components to apply continuation
! In  list_node_1    : name of list for first node of edges
! In  list_node_2    : name of list for second node of edges
! Out cmp_name       : name of component to apply continuation
! Out nb_node_sele   : final number of nodes to use in continuation
!
! --------------------------------------------------------------------------------------------------
!
    integer :: algo_lagr, i_cmp, nb_cmp
    character(len=14) :: sdline_crack
    character(len=19) :: tabai
    integer :: nb_edge
    character(len=8), pointer :: v_list_cmp(:) => null()
    aster_logical :: l_pilo, l_ainter
!
! --------------------------------------------------------------------------------------------------
!
    tabai = '&&NMMEIN.TABAI'
    l_ainter = .false.
    sdline_crack = '&&NMMEIN.LISEQ'
    algo_lagr = 2
    nb_node_sele = 0
    l_pilo = .true.
    call jelira(list_cmp, 'LONMAX', ival=nb_cmp)
!
! - Lagrange multiplier space selection for contact
!
    call xlagsp(mesh, model, crack, algo_lagr, nb_dim, &
                sdline_crack, l_pilo, tabai, l_ainter)
!
! - Init continuation method for XFEM
!
    call jelira(sdline_crack, 'LONMAX', ival=nb_edge)
    nb_edge = nb_edge/2
    call nmaret(nb_edge, nb_node, nb_dim, sdline_crack, nb_node_sele, &
                list_node, list_node_1, list_node_2)
!
! - Modification of list of components
!
    call jeveuo(list_cmp, 'E', vk8=v_list_cmp)
    do i_cmp = 1, nb_cmp
        cmp_name = v_list_cmp(i_cmp)
        if (cmp_name .eq. 'DX') v_list_cmp(i_cmp) = 'H1X'
        if (cmp_name .eq. 'DY') v_list_cmp(i_cmp) = 'H1Y'
        if (cmp_name .eq. 'DZ') v_list_cmp(i_cmp) = 'H1Z'
        if (cmp_name(1:4) .eq. 'DTAN' .or. cmp_name .eq. 'DNOR') then
            call jedetr(list_cmp)
            call wkvect(list_cmp, 'V V K8', nb_dim, vk8=v_list_cmp)
            v_list_cmp(1) = 'H1X'
            v_list_cmp(2) = 'H1Y'
            if (nb_dim .eq. 3) v_list_cmp(3) = 'H1Z'
            goto 2
        end if
    end do
2   continue
!
    call jedetr(sdline_crack)
end subroutine
