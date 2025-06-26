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

subroutine afchno(chamn, base, gran_name, mesh, nb_node, &
                  nbcpno, desc, nb_equa, typval, rval, &
                  cval, kval)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
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
    integer(kind=8) :: nbcpno(*), desc(*)
    real(kind=8), optional :: rval(*)
    complex(kind=8), optional :: cval(*)
    character(len=*), optional :: kval(*)
    character(len=*) :: chamn, gran_name, base, typval, mesh
!
!
!
    character(len=19) :: chamno, nume_equa
    integer(kind=8) :: ncmp, ncmpmx
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, ic, idec, iec, ii, inec
    integer(kind=8) :: ino, jj, lnueq, nb_equa, lvale, nb_node
    integer(kind=8) :: nec, nn, idx_gd
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
    aster_logical :: l_affe
!-----------------------------------------------------------------------
    call jemarq()
    chamno = chamn
    l_affe = present(rval) .or. present(cval) .or. present(kval)
    if (l_affe) then
        if (typval(1:1) .eq. 'R') then
            ASSERT(present(rval))
        else if (typval(1:1) .eq. 'C') then
            ASSERT(present(cval))
        else if (typval(1:2) .eq. 'K8') then
            ASSERT(present(kval))
        end if
    end if
!
    call jenonu(jexnom('&CATA.GD.NOMGD', gran_name), idx_gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', idx_gd), 'LONMAX', ncmpmx)
    call dismoi('NB_EC', gran_name, 'GRANDEUR', repi=nec)
!
! - Create NUME_EQUA
!
    nume_equa = chamno(1:8)//'.NUME_EQUA '
    call crprno(nume_equa, base, mesh, gran_name, nb_equa)
!
! - Create NODE field
!
    call vtcreb(chamno, base, typval, &
                meshz=mesh, nume_equaz=nume_equa, idx_gdz=idx_gd, nb_equa_inz=nb_equa)
!
!     --- AFFECTATION DU .PRNO DE L'OBJET NUME_EQUA ---
!
    call jeveuo(nume_equa//'.PRNO', 'E', vi=prno)
    ii = 0
    idec = 1
    do ino = 1, nb_node
        prno((nec+2)*(ino-1)+1) = idec
        prno((nec+2)*(ino-1)+2) = nbcpno(ino)
        do inec = 1, nec
            ii = ii+1
            prno((nec+2)*(ino-1)+2+inec) = desc(ii)
        end do
        idec = idec+nbcpno(ino)
    end do
!
!     --- AFFECTATION DU .VALE DE L'OBJET CHAMNO ---
!
    if (l_affe) then
        call jeveuo(chamno//'.VALE', 'E', lvale)
        call jeveuo(nume_equa//'.NUEQ', 'E', lnueq)
        do ino = 1, nb_node
            i1 = prno((nec+2)*(ino-1)+1)+lnueq-1
            do ic = 1, ncmpmx
                iec = (ic-1)/30+1
                jj = ic-30*(iec-1)
                ii = 2**jj
                nn = iand(desc((ino-1)*nec+iec), ii)
                if (nn .gt. 0) then
                    if (typval(1:1) .eq. 'R') then
                        zr(lvale-1+zi(i1)) = rval((ino-1)*ncmpmx+ic)
                    else if (typval(1:1) .eq. 'C') then
                        zc(lvale-1+zi(i1)) = cval((ino-1)*ncmpmx+ic)
                    else if (typval(1:2) .eq. 'K8') then
                        zk8(lvale-1+zi(i1)) = kval((ino-1)*ncmpmx+ic)
                    end if
                    i1 = i1+1
                end if
            end do
        end do
    end if
!
! - Create object local components (field) => global components (catalog)
!
    call cmpcha(chamno, cmp_name, cata_to_field, field_to_cata, nb_cmpz=ncmp)
!
! - Compute .DEEQ object
!
    call pteequ(nume_equa, base, nb_equa, idx_gd, ncmp, &
                field_to_cata)
    AS_DEALLOCATE(vi=cata_to_field)
    AS_DEALLOCATE(vi=field_to_cata)
    AS_DEALLOCATE(vk8=cmp_name)
!
    call jedema()
end subroutine
