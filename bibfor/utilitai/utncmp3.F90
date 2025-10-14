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

subroutine utncmp3(nume_equa, ncmp, list_cmp, list_name)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: ncmp
    character(len=*) :: nume_equa, list_cmp, list_name
!
!     RENVOIE UN LIEN iCMP -> NOM
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: jprno, gd, nec, tabec(11), j, ino, iec, icmp, ncmpmx
    integer(kind=8) ::  iad, kcmp, nnoe, ncmpt
    integer(kind=8) :: jcmp
    character(len=8) :: noma, exilag
    character(len=19) :: nume_equa19
    integer(kind=8), pointer :: vicmp(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    nume_equa19 = nume_equa
    ncmp = 0
!
    call dismoi('NOM_MAILLA', nume_equa19, 'NUME_EQUA', repk=noma)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nnoe)

    call dismoi('NUM_GD', nume_equa19, 'NUME_EQUA', repi=gd)
!
    nec = nbec(gd)
    ASSERT(nec .le. 11)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad)
    AS_ALLOCATE(vi=vicmp, size=ncmpmx)
!
!     ==================================================================
!                            NUME_EQUA
!     ==================================================================

    call jeveuo(jexnum(nume_equa19//'.PRNO', 1), 'L', jprno)
    do ino = 1, nnoe
        do iec = 1, nec
            tabec(iec) = zi(jprno-1+(ino-1)*(nec+2)+2+iec)
        end do
        do icmp = 1, ncmpmx
            if (exisdg(tabec, icmp)) then
                do j = 1, ncmp
                    if (vicmp(j) .eq. icmp) goto 14
                end do
                ncmp = ncmp+1
                vicmp(ncmp) = icmp
            end if
14          continue
        end do
    end do
!
    if (ncmp .eq. 0) then
        call utmess('F', 'UTILITAI5_53')
    end if
!
    call dismoi('EXIS_LAGR', nume_equa19, 'NUME_EQUA', repk=exilag)
    ncmpt = ncmp
    print *, exilag
    if (exilag .eq. 'OUI') then
        ncmpt = ncmpt+1
        print *, ncmp, ncmpt
    end if

!
    call wkvect(list_name, 'V V K8', ncmpt, kcmp)
    call wkvect(list_cmp, 'V V I', ncmpt, jcmp)
    do icmp = 1, ncmp
        zk8(kcmp+icmp-1) = zk8(iad-1+vicmp(icmp))
        zi(jcmp+icmp-1) = vicmp(icmp)
    end do
    if (exilag .eq. 'OUI') then
        zk8(kcmp+ncmpt-1) = "LAGR"
        zi(jcmp+ncmpt-1) = 0
    end if
    ncmp = ncmpt
    AS_DEALLOCATE(vi=vicmp)
!
    call jedema()
end subroutine
