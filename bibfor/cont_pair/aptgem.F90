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

subroutine aptgem(sdappa, mesh, newgeo, sdcont_defi, model_ndim, &
                  i_zone, zone_type, epsi_maxi, jdecma, &
                  nb_elem, err_appa)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterc/r8maem.h"
#include "asterfort/apcoma.h"
#include "asterfort/apcond.h"
#include "asterfort/apcpoi.h"
#include "asterfort/apcpou.h"
#include "asterfort/cfnben.h"
#include "asterfort/cfnumm.h"
#include "asterfort/aptypm.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/mmctan.h"
#include "asterfort/mmtann.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sdcont_defi
    character(len=19), intent(in) :: newgeo
    integer(kind=8), intent(in) :: model_ndim
    integer(kind=8), intent(in) :: i_zone
    integer(kind=8), intent(in) :: jdecma
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(inout) :: err_appa
    character(len=4), intent(in) :: zone_type
    real(kind=8), intent(in) :: epsi_maxi
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Compute tangents at each node for each element - on current zone
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  mesh             : name of mesh
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  model_ndim       : dimension of model
! In  i_zone           : index of contact zone
! In  jdecma           : shift in contact datastructure for the beginning of element in contact zone
! In  nb_elem          : number of elements in contact zone
! In  zone_type        : type of zone
!                        'MAIT' for master
!                        'ESCL' for slave
! In  epsi_maxi        : maximum tolerance for Newton algorithm
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elem_type, elem_name, node_name, valk(2)
    mpi_int :: i_proc, nb_proc, mpicou
    integer(kind=8) :: nb_elem_mpi, nbr_elem_mpi, idx_start, idx_end
    integer(kind=8) :: node_nume(9), longc
    integer(kind=8) :: elem_nbnode, niverr
    aster_logical :: l_beam, l_poi1, one_proc
    integer(kind=8) :: i_node, i_elem, elem_ndim
    integer(kind=8) :: elem_indx, elem_nume
    real(kind=8) :: tau1(3), tau2(3)
    character(len=24) :: sdappa_tgel
    real(kind=8) :: elem_coor(27), node_coor(3)
    real(kind=8), pointer :: v_sdappa_tgel(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    one_proc = .false._1
!
! - Acces to pairing datastructure
!
    sdappa_tgel = sdappa(1:19)//'.TGEL'
!
! - MPI informations
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(mpicou, rank=i_proc, size=nb_proc)
    if (one_proc) then
        nb_proc = 1
        i_proc = 0
    end if
    nb_elem_mpi = int(nb_elem/nb_proc)
    nbr_elem_mpi = nb_elem-nb_elem_mpi*nb_proc
    idx_start = 1+(i_proc)*nb_elem_mpi
    idx_end = idx_start+nb_elem_mpi-1+nbr_elem_mpi*int((i_proc+1)/nb_proc)
!
! - Loop on elements
!
    do i_elem = idx_start, idx_end
!
! ----- Current element
!
        elem_indx = i_elem+jdecma
        call cfnumm(sdcont_defi, elem_indx, elem_nume)
!
! ----- Number of nodes
!
        call cfnben(sdcont_defi, elem_indx, 'CONNEX', enti_nb_=elem_nbnode)
!
! ----- Parameters of current element
!
        call aptypm(mesh, elem_nume, elem_ndim, elem_nbnode, elem_type, &
                    elem_name)
!
! ----- Coordinates of current element
!
        call apcoma(mesh, newgeo, elem_nume, elem_nbnode, elem_coor)
!
! ----- Get absolute index of nodes
!
        call jeveuo(jexnum(mesh//'.CONNEX', elem_nume), 'L', vi=v_mesh_connex)
        do i_node = 1, elem_nbnode
            node_nume(i_node) = v_mesh_connex(i_node)
        end do
!
! ----- Right length
!
        call jelira(jexnum(sdappa_tgel, elem_indx), 'LONMAX', longc)
        longc = longc/6
!
! ----- Special types of element
!
        l_beam = (elem_type(1:2) .eq. 'SE') .and. (model_ndim .eq. 3)
        l_poi1 = elem_type .eq. 'PO1'
!
! ----- Current element
!
        call jeveuo(jexnum(sdappa_tgel, elem_indx), 'E', vr=v_sdappa_tgel)
!
! ----- Loop on nodes
!
        do i_node = 1, elem_nbnode
            tau1(1:3) = r8maem()
            tau2(1:3) = r8maem()
!
! --------- Current node
!
            call apcond(newgeo, node_nume(i_node), node_coor)
            node_name = int_to_char8(node_nume(i_node))
!
! --------- Compute local basis for these node
!
            if (l_poi1) then
                call apcpoi(sdcont_defi, model_ndim, i_zone, elem_name, &
                            zone_type, tau1, tau2)
            else
                call mmctan(elem_nume, elem_type, elem_nbnode, elem_ndim, elem_coor, &
                            node_coor, epsi_maxi, tau1, tau2, err_appa)
                if (l_beam) then
                    call apcpou(sdcont_defi, i_zone, elem_name, zone_type, &
                                tau1, tau2)
                end if
            end if
!
! --------- Norm
!
            call mmtann(model_ndim, tau1, tau2, niverr)
            if (niverr .eq. 1) then
                valk(1) = elem_name
                valk(2) = node_name
                call utmess('F', 'APPARIEMENT_14', nk=2, valk=valk)
            end if
!
! --------- Save tangents (careful: for QUAD8 only 4 nodes in contact datastructures!)
!
            if (i_node .le. longc) then
                v_sdappa_tgel(6*(i_node-1)+1) = tau1(1)
                v_sdappa_tgel(6*(i_node-1)+2) = tau1(2)
                v_sdappa_tgel(6*(i_node-1)+3) = tau1(3)
                v_sdappa_tgel(6*(i_node-1)+4) = tau2(1)
                v_sdappa_tgel(6*(i_node-1)+5) = tau2(2)
                v_sdappa_tgel(6*(i_node-1)+6) = tau2(3)
            end if
        end do
    end do
!
    call jedema()
end subroutine
