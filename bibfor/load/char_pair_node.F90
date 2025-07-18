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

subroutine char_pair_node(mesh, nb_node, &
                          list_node_i1, list_node_i2, list_node_o1, list_node_o2, i_error)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/matrot.h"
#include "asterfort/padist.h"
#include "asterfort/parotr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
! Person in charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nb_node
    character(len=24), intent(in) :: list_node_i1
    character(len=24), intent(in) :: list_node_i2
    character(len=24), intent(in) :: list_node_o1
    character(len=24), intent(in) :: list_node_o2
    integer(kind=8), intent(out) :: i_error
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Pairing two list of nodes with transformation
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh          : name of mesh
! In  nb_node       : number of nodes read
! In  list_node_i1  : first list of nodes to pair
! In  list_node_i2  : second list of nodes to pair
! In  list_node_o1  : first list of nodes to pair after paring
! In  list_node_o2  : second list of nodes to pair after paring
! Out i_error       : 0 if OK, 1 if errors occurred
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ino_1, ino_2, i_no
    integer(kind=8) :: j_node_i1, j_node_i2
    integer(kind=8) :: j_node_o1, j_node_o2, j_node_o3, j_node_o4
    integer(kind=8) :: j_node_inv
    integer(kind=8) :: jgeom_init
    integer(kind=8) :: ier, iexcor
    integer(kind=8) :: nume_node_1, nume_node_2, nume_node_3, nume_node_4, nume_node_a
    character(len=8) :: name_node_1, name_node_2, name_node_3, name_node_4, name_node_a
    integer(kind=8) :: ino_mini
    real(kind=8) :: dist, dist_mini, cent(3), tran(3)
    real(kind=8) :: matr_rota(3, 3), x1(3), x2(3)
    character(len=24) :: valk(3)
    real(kind=8) :: angl_naut_d(3)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ier = 0
    i_error = 0
    cent(1:3) = 0.d0
    tran(1:3) = 0.d0
    call jeveuo(mesh//'.COORDO    .VALE', 'L', jgeom_init)
    angl_naut_d(1:3) = 0.d0
!
! - Rotation matrix
!
    call matrot(angl_naut_d, matr_rota)
!
! - Input lists
!
    call jeveuo(list_node_i1, 'L', j_node_i1)
    call jeveuo(list_node_i2, 'L', j_node_i2)
!
! - Output lists
!
    call jeveuo(list_node_o1, 'L', j_node_o1)
    call jeveuo(list_node_o2, 'L', j_node_o2)
!
! - Working objects
!
    call wkvect('&&PACOAP.LISOU3', 'V V I', nb_node, j_node_o3)
    call wkvect('&&PACOAP.LISOU4', 'V V I', nb_node, j_node_o4)
    call wkvect('&&PACOAP.LISINV', 'V V I', nb_node, j_node_inv)
!
! - Pairing: list_node_1 -> list_node_2
!
    do ino_1 = 1, nb_node
!
        nume_node_1 = zi(j_node_i1-1+ino_1)
!
! ----- Apply transformation for translation/rotation
!
        call parotr(mesh, jgeom_init, nume_node_1, 0, cent, &
                    matr_rota, tran, x1)
!
! ----- Find nearest node
!
        dist_mini = r8gaem()
        ino_mini = 0
        do ino_2 = 1, nb_node
            nume_node_2 = zi(j_node_i2-1+ino_2)
            x2(1) = zr(jgeom_init+3*(nume_node_2-1)-1+1)
            x2(2) = zr(jgeom_init+3*(nume_node_2-1)-1+2)
            x2(3) = zr(jgeom_init+3*(nume_node_2-1)-1+3)
            dist = padist(3, x1, x2)
            if (dist .lt. dist_mini) then
                dist_mini = dist
                ino_mini = ino_2
            end if
        end do
!
! ----- No nearest node found
!
        if (ino_mini .eq. 0) then
            i_error = 1
            name_node_1 = int_to_char8(nume_node_1)
            call utmess('E', 'CHARGES2_50', sk=name_node_1)
            goto 99
        end if
        nume_node_2 = zi(j_node_i2-1+ino_mini)
!
! ----- Nearest node found
!
        nume_node_2 = zi(j_node_i2-1+ino_mini)
        nume_node_a = zi(j_node_inv-1+ino_mini)
        if (nume_node_a .eq. 0) then
            zi(j_node_o1-1+ino_1) = nume_node_1
            zi(j_node_o2-1+ino_1) = nume_node_2
            zi(j_node_inv-1+ino_mini) = nume_node_2
        else
            ier = ier+1
            name_node_1 = int_to_char8(nume_node_1)
            name_node_2 = int_to_char8(nume_node_2)
            name_node_a = int_to_char8(nume_node_a)
            valk(1) = name_node_2
            valk(2) = name_node_1
            valk(3) = name_node_a
            call utmess('E', 'CHARGES2_51', nk=3, valk=valk)
        end if
    end do
!
! - Conflict: two nodes have same neighbour
!
    if (ier .ne. 0) then
        i_error = 1
        goto 99
    end if
!
    do i_no = 1, nb_node
        zi(j_node_inv-1+i_no) = 0
    end do
!
! - Pairing: list_node_2 -> list_node_1
!
    do ino_2 = 1, nb_node
        nume_node_2 = zi(j_node_i2-1+ino_2)
!
        x2(1) = zr(jgeom_init-1+3*(nume_node_2-1)+1)
        x2(2) = zr(jgeom_init-1+3*(nume_node_2-1)+2)
        x2(3) = zr(jgeom_init-1+3*(nume_node_2-1)+3)
!
! ----- Find nearest node
!
        dist_mini = r8gaem()
        ino_mini = 0
        do ino_1 = 1, nb_node
            nume_node_1 = zi(j_node_i1-1+ino_1)
!
! --------- Apply transformation for translation/rotation
!
            call parotr(mesh, jgeom_init, nume_node_1, 0, cent, &
                        matr_rota, tran, x1)
            dist = padist(3, x1, x2)
            if (dist .lt. dist_mini) then
                dist_mini = dist
                ino_mini = ino_1
            end if
        end do
!
! ----- No nearest node found
!
        if (ino_mini .eq. 0) then
            i_error = 1
            name_node_2 = int_to_char8(nume_node_2)
            call utmess('F', 'CHARGES2_50', sk=name_node_2)
            goto 99
        end if
        nume_node_1 = zi(j_node_i1-1+ino_mini)
!
! ----- Nearest node found
!
        nume_node_1 = zi(j_node_i1-1+ino_mini)
        nume_node_a = zi(j_node_inv-1+ino_mini)
        if (nume_node_a .eq. 0) then
            zi(j_node_o3-1+ino_2) = nume_node_2
            zi(j_node_o4-1+ino_2) = nume_node_1
            zi(j_node_inv-1+ino_mini) = nume_node_1
        else
            ier = ier+1
            name_node_1 = int_to_char8(nume_node_1)
            name_node_2 = int_to_char8(nume_node_2)
            name_node_a = int_to_char8(nume_node_a)
            valk(1) = name_node_1
            valk(2) = name_node_2
            valk(3) = name_node_a
            call utmess('E', 'CHARGES2_51', nk=3, valk=valk)
        end if
    end do
!
! - Conflict: two nodes have same neighbour
!
    if (ier .ne. 0) then
        i_error = 1
        goto 99
    end if
!
! - Check: (list_node_1 -> list_node_2) == (list_node_2 -> list_node_1)
!
    do ino_1 = 1, nb_node
        iexcor = 0
        nume_node_1 = zi(j_node_o1-1+ino_1)
        nume_node_2 = zi(j_node_o2-1+ino_1)
        name_node_1 = int_to_char8(nume_node_1)
        name_node_2 = int_to_char8(nume_node_2)
        do ino_2 = 1, nb_node
            nume_node_3 = zi(j_node_o3-1+ino_2)
            nume_node_4 = zi(j_node_o4-1+ino_2)
            name_node_3 = int_to_char8(nume_node_3)
            name_node_4 = int_to_char8(nume_node_4)
            if (nume_node_2 .eq. nume_node_3) then
                iexcor = 1
                if (nume_node_1 .ne. nume_node_4) then
                    ier = ier+1
                    valk(1) = name_node_1
                    valk(2) = name_node_2
                    valk(3) = name_node_4
                    call utmess('E', 'CHARGES2_51', nk=3, valk=valk)
                end if
            end if
        end do
        if (iexcor .eq. 0) then
            ier = ier+1
            valk(1) = name_node_1
            call utmess('E', 'CHARGES2_52', sk=valk(1))
        end if
    end do
!
    if (ier .ne. 0) then
        i_error = 1
        goto 99
    end if
!
99  continue
!
    call jedetr('&&PACOAP.LISOU3')
    call jedetr('&&PACOAP.LISOU4')
    call jedetr('&&PACOAP.LISINV')
!
    call jedema()
end subroutine
