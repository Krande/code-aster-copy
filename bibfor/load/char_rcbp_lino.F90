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

subroutine char_rcbp_lino(mesh, name_ancr, list_node, nb_node)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jelira.h"
#include "asterfort/wkvect.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: name_ancr
    character(len=24), intent(in) :: list_node
    integer(kind=8), intent(out) :: nb_node
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! RELA_CINE_BP - Get list of nodes for ancrage
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh         : name of mesh
! In  name_ancr    : name of ancrage
! In  list_node    : list of nodes of ancrage
! Out nb_node      : number of nodes of ancrage
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jgro, jlino_old, jind, jlino_new
    character(len=8) :: k8bid
    integer(kind=8) :: ino, nbno, numnoe, indnoe, ino_1, indlis
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    nb_node = 0
!
! - Acces to nodes
!
    call jeveuo(jexnom(mesh//'.GROUPENO', name_ancr), 'L', jgro)
    call jelira(jexnom(mesh//'.GROUPENO', name_ancr), 'LONUTI', nbno, k8bid)
!
! - Create list of nodes
!
    call wkvect('&&CAPREC.LIST', 'V V I', nbno, jlino_old)
    indnoe = 0
    do ino = 1, nbno
        numnoe = zi(jgro+ino-1)
        indnoe = indnoe+1
        zi(jlino_old+indnoe-1) = numnoe
    end do
!
! - No double
!
    call wkvect('&&CAPREC.INDICE', 'V V I', nbno, jind)
    do ino = 1, nbno
        do ino_1 = ino+1, nbno
            if (zi(jlino_old+ino_1-1) .eq. zi(jlino_old+ino-1)) then
                zi(jind+ino_1-1) = 1
            end if
        end do
    end do
    call wkvect(list_node, 'V V I', nbno, jlino_new)
!
! - Re-create list of nodes
!
    indlis = 0
    do ino = 1, nbno
        if (zi(jind+ino-1) .eq. 0) then
            indlis = indlis+1
            zi(jlino_new+indlis-1) = zi(jlino_old+ino-1)
        end if
    end do
    nb_node = indlis
!
    call jedetr('&&CAPREC.INDICE')
    call jedetr('&&CAPREC.LIST')
!
    call jedema()
end subroutine
