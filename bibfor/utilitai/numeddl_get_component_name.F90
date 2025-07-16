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

subroutine numeddl_get_component_name(nume19, cmpid, cmpname)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: cmpid
    character(len=*) :: nume19
    character(len=8) :: cmpname
!
!     Returns a list of all available components (aka dof) in the
!     numbering or for a given node
!     -----------------------------------------------------------------------------
!     In    nume19      : name of numbering
!     In    allorone    : "ONE" or "ALL" to get components for one or all nodes
!     In    nodeid      : id of the node (1-based number)
!     Out   ncmp        : number of components
!     Out   stringarray : name  of components
!     In    maxcmp      : max possible number of components (for safety purpose)
!     -----------------------------------------------------------------------------
!
    integer(kind=8) :: nb_cmp_gd, vali(2), jcmp
    character(len=8) :: nomgd
    character(len=19) :: numeddl
    character(len=24), pointer :: refn(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    numeddl = nume19
    call jeveuo(numeddl(1:14)//'.NUME.REFN', 'L', vk24=refn)
    nomgd = refn(2) (1:8)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmp)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', nb_cmp_gd)
    if ((cmpid .gt. nb_cmp_gd) .or. (cmpid .le. 0)) then
        vali(1) = cmpid
        vali(2) = nb_cmp_gd
        call utmess('F', 'UTILITAI_30', ni=2, vali=vali)
    end if
    cmpname = zk8(jcmp-1+cmpid)
!
    call jedema()
end subroutine
