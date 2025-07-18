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
subroutine apverl(sdappa, mesh, sdcont_defi)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterc/r8rddg.h"
#include "asterfort/cfinvm.h"
#include "asterfort/cfnben.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfnumn.h"
#include "asterfort/cfdisi.h"
#include "asterfort/mminfi.h"
#include "asterfort/cfcald.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mmnorm.h"
#include "blas/ddot.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sdcont_defi
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Check normals discontinuity
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  mesh             : name of mesh
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: node_mast_name, elem_mast_name, node_old
    integer(kind=8) :: nb_cont_zone, model_ndim
    integer(kind=8) :: i_zone
    integer(kind=8) :: jdecnm, nb_node_mast
    integer(kind=8) :: i_node_mast, i_elem, i_node_curr, i_elem_node
    integer(kind=8) :: jdeciv
    integer(kind=8) :: elem_mast_indx, elem_mast_nume, node_mast_indx(1), node_mast_nume(1)
    integer(kind=8) :: node_nbelem, elem_nbnode
    real(kind=8) :: tau1(3), tau2(3), norm(3)
    real(kind=8) :: tau1_node(3), tau2_node(3), normnd(3)
    real(kind=8) :: noor1, noor2
    real(kind=8) :: angmax
    character(len=24) :: sdappa_tgel, sdappa_tgno
    real(kind=8), pointer :: v_sdappa_tgel(:) => null()
    real(kind=8), pointer :: v_sdappa_tgno(:) => null()
    character(len=24) :: sdappa_verk
    character(len=8), pointer :: v_sdappa_verk(:) => null()
    character(len=24) :: sdappa_vera
    real(kind=8), pointer :: v_sdappa_vera(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    real(kind=8) :: prosca, angle, angl_old, val
    integer(kind=8) :: inoeu
    aster_logical :: apcald
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Get contact datastructures
!
    sdappa_verk = sdappa(1:19)//'.VERK'
    sdappa_vera = sdappa(1:19)//'.VERA'
    call jeveuo(sdappa_verk, 'E', vk8=v_sdappa_verk)
    call jeveuo(sdappa_vera, 'E', vr=v_sdappa_vera)
    call jelira(sdappa_verk, 'LONUTI', inoeu)
    if (inoeu .ne. 0) goto 999
    angmax = 5.d0
!
! - Acces to pairing datastructure
!
    sdappa_tgel = sdappa(1:19)//'.TGEL'
    sdappa_tgno = sdappa(1:19)//'.TGNO'
    call jeveuo(sdappa_tgno, 'L', vr=v_sdappa_tgno)
!
! - Get parameters
!
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
!
! - Loop on contact zones
!
    inoeu = 0
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters on current zone - Master
!
        nb_node_mast = mminfi(sdcont_defi, 'NBNOM', i_zone)
        jdecnm = mminfi(sdcont_defi, 'JDECNM', i_zone)
        apcald = cfcald(sdcont_defi, i_zone, 'MAIT')
        if (apcald) then
!
! --------- Loop on nodes
!
            do i_node_mast = 1, nb_node_mast
!
! ------------- Current node
!
                node_mast_indx(1) = i_node_mast+jdecnm
                call cfnumn(sdcont_defi, 1, node_mast_indx(1), node_mast_nume(1))
                node_mast_name = int_to_char8(node_mast_nume(1))
!
! ------------- Get tangents
!
                tau1_node(1) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+1)
                tau1_node(2) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+2)
                tau1_node(3) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+3)
                tau2_node(1) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+4)
                tau2_node(2) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+5)
                tau2_node(3) = v_sdappa_tgno(6*(node_mast_indx(1)-1)+6)
!
! ------------- Compute normal
!
                call mmnorm(model_ndim, tau1_node, tau2_node, normnd, noor2)
!
! ------------- Number of elements attached to node
!
                call cfnben(sdcont_defi, node_mast_indx(1), 'CONINV', node_nbelem, jdeciv)
!
! ------------- Loop on elements attached to node
!
                do i_elem = 1, node_nbelem
!
! ----------------- Current element
!
                    call cfinvm(sdcont_defi, jdeciv, i_elem, elem_mast_indx)
                    call cfnumm(sdcont_defi, elem_mast_indx, elem_mast_nume)
                    elem_mast_name = int_to_char8(elem_mast_nume)
                    call cfnben(sdcont_defi, elem_mast_indx, 'CONNEX', elem_nbnode)
!
! ----------------- Access to connectivity
!
                    call jeveuo(jexnum(mesh//'.CONNEX', elem_mast_nume), 'L', vi=v_mesh_connex)
!
! ----------------- Get current index of node
!
                    i_node_curr = 0
                    do i_elem_node = 1, elem_nbnode
                        if (v_mesh_connex(i_elem_node) .eq. node_mast_nume(1)) then
                            i_node_curr = i_elem_node
                        end if
                    end do
                    ASSERT(i_node_curr .ne. 0)
!
! ----------------- Access to current tangent
!
                    call jeveuo(jexnum(sdappa_tgel, elem_mast_indx), 'L', vr=v_sdappa_tgel)
                    tau1(1) = v_sdappa_tgel(6*(i_node_curr-1)+1)
                    tau1(2) = v_sdappa_tgel(6*(i_node_curr-1)+2)
                    tau1(3) = v_sdappa_tgel(6*(i_node_curr-1)+3)
                    tau2(1) = v_sdappa_tgel(6*(i_node_curr-1)+4)
                    tau2(2) = v_sdappa_tgel(6*(i_node_curr-1)+5)
                    tau2(3) = v_sdappa_tgel(6*(i_node_curr-1)+6)
!
! ----------------- Compute normal
!
                    call mmnorm(model_ndim, tau1, tau2, norm, noor1)
!
! ----------------- Compute angle
!
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    prosca = ddot(b_n, norm, b_incx, normnd, b_incy)
                    if (abs(noor1*noor2) .gt. r8prem()) then
                        val = prosca/(noor1*noor2)
                        if (val .gt. 1.d0) then
                            val = 1.d0
                        end if
                        if (val .lt. -1.d0) then
                            val = -1.d0
                        end if
                        angle = acos(val)
                        angle = angle*r8rddg()
                        angl_old = v_sdappa_vera(node_mast_nume(1))
                        node_old = v_sdappa_verk(node_mast_nume(1))
                        if (angle .gt. angmax) then
                            if (node_old .eq. ' ') then
                                inoeu = inoeu+1
                                if (angl_old .lt. angle) then
                                    v_sdappa_verk(node_mast_nume(1)) = node_mast_name
                                    v_sdappa_vera(node_mast_nume(1)) = angle
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end if
    end do
!
    call jeecra(sdappa_verk, 'LONUTI', inoeu)
!
999 continue
!
    call jedema()
!
end subroutine
