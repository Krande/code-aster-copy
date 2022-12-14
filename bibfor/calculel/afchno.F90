! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine afchno(chamn, base, gran_name, mesh, nb_node,&
                  nbcpno, desc, nb_equa, typval, rval,&
                  cval, kval)
    implicit none
#include "jeveux.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cmpcha.h"
#include "asterfort/vtcreb.h"
#include "asterfort/crprno.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/pteequ.h"
!
    integer :: nbcpno(*), desc(*)
    real(kind=8) :: rval(*)
    complex(kind=8) :: cval(*)
    character(len=*) :: chamn, gran_name, base, typval, kval(*), mesh
!
!
!
    character(len=19) :: chamno, prof_chno
    integer :: ncmp, ncmpmx
!
!-----------------------------------------------------------------------
    integer :: i1, ic, idec, iec, ii, inec
    integer :: ino, jj, lnueq, nb_equa, lvale, nb_node
    integer :: nec, nn, idx_gd
    integer, pointer :: cata_to_field(:) => null()
    integer, pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
    integer, pointer :: prno(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    chamno = chamn
!
    call jenonu(jexnom('&CATA.GD.NOMGD', gran_name), idx_gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', idx_gd), 'LONMAX', ncmpmx)
    call dismoi('NB_EC', gran_name, 'GRANDEUR', repi=nec)
!
! - Create PROF_CHNO
!
    prof_chno = chamno(1:8)//'.PROF_CHNO '
    call crprno(prof_chno, base, mesh, gran_name, nb_equa)
!
! - Create NODE field
!
    call vtcreb(chamno      , base                  , typval,&
                meshz = mesh, prof_chnoz = prof_chno, idx_gdz = idx_gd, nb_equa_inz = nb_equa)
!
!     --- AFFECTATION DU .PRNO DE L'OBJET PROF_CHNO ---
!
    call jeveuo(prof_chno//'.PRNO', 'E', vi=prno)
    ii = 0
    idec = 1
    do ino = 1, nb_node
        prno((nec+2)*(ino-1)+1) = idec
        prno((nec+2)*(ino-1)+2) = nbcpno(ino)
        do inec = 1, nec
            ii = ii + 1
            prno((nec+2)*(ino-1)+2+inec) = desc(ii)
        end do
        idec = idec + nbcpno(ino)
    end do
!
!     --- AFFECTATION DU .VALE DE L'OBJET CHAMNO ---
!
    call jeveuo(chamno//'.VALE', 'E', lvale)
    call jeveuo(prof_chno//'.NUEQ', 'E', lnueq)
    do ino = 1, nb_node
        i1 = prno((nec+2)*(ino-1)+1) + lnueq - 1
        do ic = 1, ncmpmx
            iec = ( ic - 1 ) / 30 + 1
            jj = ic - 30 * ( iec - 1 )
            ii = 2**jj
            nn = iand( desc((ino-1)*nec+iec) , ii )
            if (nn .gt. 0) then
                if (typval(1:1) .eq. 'R') then
                    zr(lvale-1+zi(i1)) = rval((ino-1)*ncmpmx+ic)
                else if (typval(1:1).eq.'C') then
                    zc(lvale-1+zi(i1)) = cval((ino-1)*ncmpmx+ic)
                else if (typval(1:2).eq.'K8') then
                    zk8(lvale-1+zi(i1)) = kval((ino-1)*ncmpmx+ic)
                endif
                i1 = i1 + 1
            endif
        end do
    end do
!
! - Create object local components (field) => global components (catalog)
!
    call cmpcha(chamno, cmp_name, cata_to_field, field_to_cata, nb_cmpz = ncmp)
!
! - Compute .DEEQ object
!
    call pteequ(prof_chno    , base, nb_equa, idx_gd, ncmp,&
                field_to_cata)
    AS_DEALLOCATE(vi = cata_to_field)
    AS_DEALLOCATE(vi = field_to_cata)
    AS_DEALLOCATE(vk8 = cmp_name)
!
    call jedema()
end subroutine
