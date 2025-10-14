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

subroutine numeddl_get_components(nume19, allorone, nodeid, ncmp, stringarray, maxcmp)
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
!
    integer(kind=8) :: ncmp, maxcmp, nodeid
    character(len=*) :: nume19, allorone
    character(len=8) :: stringarray(maxcmp)
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
    integer(kind=8) :: jprno, gd, nec, tabec(11), j, ino, iec, icmp, ncmpmx
    integer(kind=8) ::  iad, nnoe, idescr
    character(len=8) :: noma, exilag
    character(len=19) :: numeddl, prno
    integer(kind=8), pointer :: vicmp(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    numeddl = nume19
    ncmp = 0
!
    call dismoi('NOM_MAILLA', numeddl, 'NUME_DDL', repk=noma)
    call dismoi('NUM_GD_SI', numeddl, 'NUME_DDL', repi=gd)
!
    call jeveuo('&CATA.GD.DESCRIGD', 'L', idescr)
    nec = nbec(gd)
    ASSERT(nec .le. 11)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad)
    AS_ALLOCATE(vi=vicmp, size=ncmpmx)
!
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nnoe)
    call dismoi('NUME_EQUA', numeddl, 'NUME_DDL', repk=prno)
    ASSERT(prno .ne. ' ')
    call jeveuo(jexnum(prno//'.PRNO', 1), 'L', jprno)
    if (allorone(1:3) .eq. "ALL") then
        do ino = 1, nnoe
            do iec = 1, nec
                tabec(iec) = zi(jprno-1+(ino-1)*(nec+2)+2+iec)
            end do
            do icmp = 1, ncmpmx
                if (exisdg(tabec, icmp)) then
                    do j = 1, ncmp
                        if (vicmp(j) .eq. icmp) goto 10
                    end do
                    ncmp = ncmp+1
                    vicmp(ncmp) = icmp
                end if
10              continue
            end do
        end do
    else if (allorone(1:3) .eq. "ONE") then
        do iec = 1, nec
            tabec(iec) = zi(jprno-1+(nodeid-1)*(nec+2)+2+iec)
        end do
        do icmp = 1, ncmpmx
            if (exisdg(tabec, icmp)) then
                do j = 1, ncmp
                    if (vicmp(j) .eq. icmp) goto 20
                end do
                ncmp = ncmp+1
                vicmp(ncmp) = icmp
            end if
20          continue
        end do
    else
        ASSERT(.false.)
    end if

    if (ncmp .eq. 0) then
        call utmess('F', 'UTILITAI5_53')
    end if
    ASSERT(ncmp .le. maxcmp)
!
    do icmp = 1, ncmp
        stringarray(icmp) = zk8(iad-1+vicmp(icmp))
    end do
!   handle the case of Lagrange multipliers
    if (allorone(1:3) .eq. "ALL") then
        call dismoi('EXIS_LAGR', numeddl, 'NUME_DDL', repk=exilag)
        if (exilag .eq. 'OUI') then
            ncmp = ncmp+1
            stringarray(ncmp) = 'LAGR'
        end if
    end if
!
    AS_DEALLOCATE(vi=vicmp)
!
    call jedema()
end subroutine
