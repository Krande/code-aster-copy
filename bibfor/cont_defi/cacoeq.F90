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
subroutine cacoeq(sdcont, mesh, lLineRela, listRela)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/afrela.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmex.h"
#include "asterfort/cfsuex.h"
#include "asterfort/cncinv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: sdcont, mesh
    aster_logical, intent(out) :: lLineRela
    character(len=19), intent(out) :: listRela
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Discrete method - Create QUAD8 linear relations
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  mesh             : name of mesh
! Out lLineRela        : flag for linear relations
! Out listRela         : name of object for linear relations
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: node_nume(3), iret
    complex(kind=8) :: vale_cplx, coef_cplx(3)
    character(len=4) :: vale_type, type_coef
    character(len=8) :: vale_func, dof_name(3), node_name(3), type_name, type_name_c
    real(kind=8) :: coef_real(3), repe_defi(6), vale_real
    integer(kind=8) :: i_excl, i_zone, i_node_quad, i_elem, i_elem_c
    integer(kind=8) :: elem_nume, type_nume, repe_type(3)
    integer(kind=8) :: suppo1, suppo2, suppo3
    integer(kind=8) :: nb_cont_zone, nb_cont_node, nb_node_quad, nb_node_elem, nb_excl
    integer(kind=8), pointer :: v_list_excl(:) => null()
    character(len=19) :: connex_inv
    integer(kind=8), pointer :: v_coninv(:) => null()
    integer(kind=8), pointer :: v_coninv_longcum(:) => null()
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_noeuqu
    integer(kind=8), pointer :: v_sdcont_noeuqu(:) => null()
!
    data repe_defi/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/
!
! --------------------------------------------------------------------------------------------------
!
    vale_type = 'REEL'
    vale_func = '&FOZERO'
    vale_cplx = (0.0d0, 0.0d0)
    vale_real = 0.0d0
    coef_cplx(1) = (1.0d0, 0.0d0)
    coef_cplx(2) = (-0.5d0, 0.0d0)
    coef_cplx(2) = (-0.5d0, 0.0d0)
    coef_real(1) = 1.0d0
    coef_real(2) = -0.5d0
    coef_real(3) = -0.5d0
    repe_type(1) = 0
    repe_type(2) = 0
    repe_type(3) = 0
    type_coef = 'REEL'
    lLineRela = ASTER_FALSE
    listRela = '&&CACOEQ.RLLISTE'
    connex_inv = '&&CACOEQ.CONINV'
    nb_node_quad = 0
!
! - Access to datastructure
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_noeuqu = sdcont_defi(1:16)//'.NOEUQU'
    call jeexin(sdcont_noeuqu, iret)
    if (iret .eq. 0) then
        nb_node_quad = 0
    else
        call jeveuo(sdcont_noeuqu, 'L', vi=v_sdcont_noeuqu)
        call jelira(sdcont_noeuqu, 'LONUTI', nb_node_quad)
    end if
!
! - No QUAD8 nodes -> exit
!
    if ((iret .eq. 0) .or. (nb_node_quad .eq. 0)) then
        goto 99
    end if
!
! - Parameters
!
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    nb_cont_node = cfdisi(sdcont_defi, 'NNOCO')
!
! - Construct inverse connectivity
!
    call cncinv(mesh, [0], 0, 'V', connex_inv)
!
! - Access to mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=v_mesh_typmail)
!
! - List of nodes to suppress
!
    AS_ALLOCATE(vi=v_list_excl, size=nb_cont_node)
!
! - Loop on nodes to link
!
    nb_node_quad = nb_node_quad/3
    i_excl = 1
    nb_excl = 0
    do i_node_quad = 1, nb_node_quad
!
! ----- Get the two nodes to link with the middle one
!
        node_nume(1) = v_sdcont_noeuqu(3*(i_node_quad-1)+1)
        node_nume(2) = v_sdcont_noeuqu(3*(i_node_quad-1)+2)
        node_nume(3) = v_sdcont_noeuqu(3*(i_node_quad-1)+3)
!
! ----- Nodes to link ?
!
        if (node_nume(2) .eq. 0) then
            lLineRela = .false.
            goto 30
        else
            lLineRela = .true.
        end if
!
! ----- Name of the three nodes
!
        node_name(1) = int_to_char8(node_nume(1))
        node_name(2) = int_to_char8(node_nume(2))
        node_name(3) = int_to_char8(node_nume(3))
!
! ----- Is node to link with the middle one belongs to SANS_GROUP_NO ?
!
        do i_zone = 1, nb_cont_zone
            call cfmmex(sdcont_defi, 'CONT', i_zone, node_nume(1), suppo1)
            call cfmmex(sdcont_defi, 'CONT', i_zone, node_nume(2), suppo2)
            call cfmmex(sdcont_defi, 'CONT', i_zone, node_nume(3), suppo3)
            if ((suppo1 .eq. 1) .or. (suppo2 .eq. 1) .or. (suppo3 .eq. 1)) then
                lLineRela = .false.
                goto 30
            end if
        end do
!
! ----- Loop on elements connected to middle node
!
        call jeveuo(jexatr(connex_inv, 'LONCUM'), 'L', vi=v_coninv_longcum)
        nb_node_elem = v_coninv_longcum(node_nume(1)+1)-v_coninv_longcum(node_nume(1))
        call jeveuo(jexnum(connex_inv, node_nume(1)), 'L', vi=v_coninv)
        do i_elem = 1, nb_node_elem
!
! --------- Type of element
!
            elem_nume = v_coninv(i_elem)
            type_nume = v_mesh_typmail(elem_nume)
            call jenuno(jexnum('&CATA.TM.NOMTM', type_nume), type_name)
!
! --------- If middle node belongs to QUAD8 and TRIA*/QUAD9 -> save it to suppress in contact
!
            if (type_name(1:5) .eq. 'QUAD8') then
                do i_elem_c = 1, nb_node_elem
                    if (i_elem .ne. i_elem_c) then
                        elem_nume = v_coninv(i_elem_c)
                        type_nume = v_mesh_typmail(elem_nume)
                        call jenuno(jexnum('&CATA.TM.NOMTM', type_nume), type_name_c)
                        if ((type_name_c(1:5) .eq. 'TRIA6') .or. &
                            (type_name_c(1:5) .eq. 'TRIA7') .or. &
                            (type_name_c(1:5) .eq. 'QUAD9')) then
                            v_list_excl(i_excl) = node_nume(1)
                            i_excl = i_excl+1
                            nb_excl = nb_excl+1
                            goto 25
                        end if
                    end if
                end do
            end if
25          continue
        end do
!
! ----- Linear relations
!
        dof_name(1) = 'DX'
        dof_name(2) = 'DX'
        dof_name(3) = 'DX'
        call afrela(coef_real, coef_cplx, dof_name, node_name, repe_type, &
                    repe_defi, 3, vale_real, vale_cplx, vale_func, &
                    type_coef, vale_type, 0.d0, listRela)
        dof_name(1) = 'DY'
        dof_name(2) = 'DY'
        dof_name(3) = 'DY'
        call afrela(coef_real, coef_cplx, dof_name, node_name, repe_type, &
                    repe_defi, 3, vale_real, vale_cplx, vale_func, &
                    type_coef, vale_type, 0.d0, listRela)
        dof_name(1) = 'DZ'
        dof_name(2) = 'DZ'
        dof_name(3) = 'DZ'
        call afrela(coef_real, coef_cplx, dof_name, node_name, repe_type, &
                    repe_defi, 3, vale_real, vale_cplx, vale_func, &
                    type_coef, vale_type, 0.d0, listRela)
30      continue
    end do
!
    if (lLineRela) then
        call cfsuex(sdcont_defi, v_list_excl, nb_excl, nb_cont_zone)
        call jedetr(connex_inv)
    end if
    AS_DEALLOCATE(vi=v_list_excl)
!
99  continue
!
end subroutine
