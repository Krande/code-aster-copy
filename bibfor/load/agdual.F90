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

subroutine agdual(load, type_liai, nbrela, nbterm)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/juveca.h"
#include "asterfort/jeexin.h"
!
! person_in_load: jacques.pellet at edf.fr
!
    character(len=8), intent(in) :: load
    character(len=*), intent(in) :: type_liai
    integer(kind=8), intent(in) :: nbrela
    integer(kind=8), intent(in) :: nbterm
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Allocate/Extend char_dual datastructure for a load
!
! --------------------------------------------------------------------------------------------------
!
! In  load        : name of load
! In  type_liai   : "type" of relations : 'LIN' / 'NLIN' / '3D3' / '3D2' / ...
! In  nbrela      : number of relations of the "pack"
! In  nbterm      : number of terms of the "pack"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), pointer :: prdk(:) => null()
    character(len=8), pointer :: prdso(:) => null()
    integer(kind=8), pointer :: prdi(:) => null()
    integer(kind=8) :: nbpaqu_max, n1, n2, iexi, kpaqu, term1
    character(len=13) :: chdual
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    ASSERT(nbrela .gt. 0)
    ASSERT(nbterm .gt. 0)
    chdual = load//'.DUAL'

!   -- faut-il creer la SD ou l'agrandir ?
!   ---------------------------------------
    call jeexin(chdual//'.PRDK', iexi)
    if (iexi .eq. 0) then
        nbpaqu_max = 30
        call wkvect(load//'.DUAL.PRDK', 'G V K8', nbpaqu_max, vk8=prdk)
        call wkvect(load//'.DUAL.PRDI', 'G V I', 3*nbpaqu_max, vi=prdi)
        call wkvect(load//'.DUAL.PRDSO', 'G V K8', 4*nbpaqu_max, vk8=prdso)
        kpaqu = 1
        term1 = 1
        goto 998
    end if

!   -- faut-il agrandir la SD ?
!   ----------------------------
    call jelira(chdual//'.PRDK', 'LONMAX', ival=n1)
    call jelira(chdual//'.PRDK', 'LONUTI', ival=n2)
    if (n2 .eq. n1) then
        nbpaqu_max = 2*n1
        call juveca(load//'.DUAL.PRDK', nbpaqu_max)
        call juveca(load//'.DUAL.PRDI', 3*nbpaqu_max)
        call juveca(load//'.DUAL.PRDSO', 4*nbpaqu_max)
    end if
    call jeveuo(chdual//'.PRDK', 'E', vk8=prdk)
    call jeveuo(chdual//'.PRDI', 'E', vi=prdi)
    kpaqu = n2+1
    ASSERT(kpaqu .ge. 2)
    term1 = prdi(3*(kpaqu-2)+3)+prdi(3*(kpaqu-2)+2)

!   -- remplissage de .PRDK et .PRDI :
!   -----------------------------------
998 continue
    call jeecra(chdual//'.PRDK', 'LONUTI', ival=kpaqu)
    prdk(kpaqu) = type_liai
    prdi(3*(kpaqu-1)+1) = nbrela
    prdi(3*(kpaqu-1)+2) = nbterm
    prdi(3*(kpaqu-1)+3) = term1

    call jedema()

end subroutine
