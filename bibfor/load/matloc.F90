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

subroutine matloc(mesh, connex_inv, keywordfact, iocc, node_nume, &
                  node_name, nb_repe_elem, list_repe_elem, matr_glob_loca)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/angvx.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: connex_inv
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    character(len=8), intent(in) :: node_name
    integer(kind=8), intent(in) :: node_nume
    integer(kind=8), intent(in) :: nb_repe_elem
    integer(kind=8), intent(in) :: list_repe_elem(*)
    real(kind=8), intent(out) :: matr_glob_loca(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation - DDL_POUTRE
!
! Computation of matrix of transition from global to local coordinate system at node
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh           : meshing
! In  connex_inv     : inverse connectivity of mesh (nodes -> elements)
! In  keywordfact    : factor keyword to read
! In  iocc           : factor keyword index in AFFE_CHAR_MECA
! In  node_nume      : number of node if mesh
! In  node_name      : name of node
! In  nb_repe_elem   : number of elements of local coordinate system
! In  list_elem      : list of elements of local coordinate sysem
! Out matr_glob_loca : matrix of transition from global to local coordinate system
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nocc, i, j, jadrm, jadrvlc
    integer(kind=8) :: nb_conn_elem, elem_nume
    integer(kind=8) :: node_nume_1, node_nume_2
    real(kind=8) :: vx(3), vy(3), vz(3), vecty(3), vxn, vyn, vyp, dgrd, angl_naut(3)
    real(kind=8) :: alpha, beta, gamma
    character(len=8) :: k8bid, type_elem, elem_name, valk(2)
    character(len=24) :: connex
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    dgrd = r8dgrd()
!
! - Elements connected to node
!
    call jelira(jexnum(connex_inv, node_nume), 'LONMAX', nb_conn_elem, k8bid)
    call jeveuo(jexnum(connex_inv, node_nume), 'L', jadrm)
!
! - Access to mesh
!
    connex = mesh//'.CONNEX'
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=vale)
!
! - Get element connected to node
!
    if (nb_conn_elem .eq. 1) then
        elem_nume = zi(jadrm)
        elem_name = int_to_char8(elem_nume)
    else
        if (nb_repe_elem .eq. 0) then
            call utmess('F', 'CHARGES2_37', sk=node_name)
        end if
        if (nb_repe_elem .ne. 1) then
            call utmess('F', 'CHARGES2_38', sk=node_name)
        end if
        do i = 1, nb_repe_elem
            elem_nume = list_repe_elem(i)
            do j = 1, nb_conn_elem
                if (zi(jadrm+j-1) .eq. elem_nume) goto 24
            end do
        end do
        call utmess('F', 'CHARGES2_39', sk=node_name)
24      continue
    end if
!
! - Check element type
!
    call jenuno(jexnum('&CATA.TM.NOMTM', typmail(elem_nume)), type_elem)
    elem_name = int_to_char8(elem_nume)
    valk(2) = node_name
    valk(1) = elem_name
    if (type_elem(1:3) .ne. 'SEG') then
        call utmess('F', 'CHARGES2_40', nk=2, valk=valk)
    end if
!
! - Get nodes of element connected
!
    call jeveuo(jexnum(connex, elem_nume), 'L', jadrvlc)
    if (node_nume .eq. zi(jadrvlc)) then
        node_nume_1 = zi(jadrvlc+1)
        node_nume_2 = zi(jadrvlc)
    else
        node_nume_1 = zi(jadrvlc)
        node_nume_2 = zi(jadrvlc+1)
    end if
!
! - Construct colinear vector to element
!
    vx(1) = vale(1+3*(node_nume_2-1))-vale(1+3*(node_nume_1-1))
    vx(2) = vale(1+3*(node_nume_2-1)+1)-vale(1+3*(node_nume_1-1)+1)
    vx(3) = vale(1+3*(node_nume_2-1)+2)-vale(1+3*(node_nume_1-1)+2)
    vxn = sqrt(vx(1)**2+vx(2)**2+vx(3)**2)
    valk(2) = node_name
    valk(1) = elem_name
    if (vxn .le. r8prem()) then
        call utmess('F', 'CHARGES2_41', nk=2, valk=valk)
    end if
    vx(1) = vx(1)/vxn
    vx(2) = vx(2)/vxn
    vx(3) = vx(3)/vxn
!
! - Is VECT_Y ?
!
    call getvr8(keywordfact, 'VECT_Y', iocc=iocc, nbval=3, vect=vecty, &
                nbret=nocc)
    if (nocc .ne. 0) then
!
! ----- VECT_Y
!
        vy(1) = vecty(1)
        vy(2) = vecty(2)
        vy(3) = vecty(3)
!
! ----- Projection
!
        vyp = vx(1)*vy(1)+vx(2)*vy(2)+vx(3)*vy(3)
        vy(1) = vy(1)-vyp*vx(1)
        vy(2) = vy(2)-vyp*vx(2)
        vy(3) = vy(3)-vyp*vx(3)
        vyn = sqrt(vy(1)**2+vy(2)**2+vy(3)**2)
        vy(1) = vy(1)/vyn
        vy(2) = vy(2)/vyn
        vy(3) = vy(3)/vyn
!
! ----- Tangent vector
!
        vz(1) = vx(2)*vy(3)-vy(2)*vx(3)
        vz(2) = vx(3)*vy(1)-vy(3)*vx(1)
        vz(3) = vx(1)*vy(2)-vy(1)*vx(2)
!
! ----- Transformation matrix
!
        do i = 1, 3
            matr_glob_loca(1, i) = vx(i)
            matr_glob_loca(2, i) = vy(i)
            matr_glob_loca(3, i) = vz(i)
        end do
        goto 999
    end if
!
! --- SI ANGL_VRIL
!
    call getvr8(keywordfact, 'ANGL_VRIL', iocc=iocc, scal=gamma, nbret=nocc)
    if (nocc .ne. 0) then
!
! ----- Get nautic angles
!
        call angvx(vx, alpha, beta)
        angl_naut(1) = alpha
        angl_naut(2) = beta
        angl_naut(3) = gamma*dgrd
!
! ----- Transformation matrix
!
        call matrot(angl_naut, matr_glob_loca)
        goto 999
    end if
!
999 continue
!
    call jedema()
!
end subroutine
