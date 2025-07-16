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

subroutine nuacno(nuagez, list_nodez, chnoz)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/nueq_chck.h"
!
!
    character(len=*), intent(in) :: chnoz
    character(len=*), intent(in) :: list_nodez
    character(len=*), intent(in) :: nuagez
!
! --------------------------------------------------------------------------------------------------
!
! Convert NUAGE to CHAM_NO
!
! --------------------------------------------------------------------------------------------------
!
! In  nuagez       : name of NUAGE datastructure
! In  chnoz        : name of CHAM_NO (nodal field) datastructure
! In  list_node   : list of nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: idx_gd
    character(len=4) :: type_scal
    character(len=8) :: mesh, gran_name
    character(len=19) :: chno, list_node, nuage, numeequa
    integer(kind=8) :: iaprno, icmp
    integer(kind=8) :: icompt, i_ec, ieq, nume_pt, itype, ival, i_pt
    integer(kind=8) :: j, jnuai, jnuav, k, i_ligr_mesh
    integer(kind=8) :: kcomp, kvale, nc, ncmp, ncmpmx, nb_ec, nb_point
    integer(kind=8), pointer :: ent_cod(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: p_list_node(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    chno = chnoz
    list_node = list_nodez
    nuage = nuagez
!
    call dismoi("NUM_GD", chno, "CHAM_NO", repi=idx_gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', idx_gd), 'LONMAX', ncmpmx)
    call jenuno(jexnum('&CATA.GD.NOMGD', idx_gd), gran_name)
    nb_ec = nbec(idx_gd)
    call wkvect('&&NUACNO.NOMCMP', 'V V I', ncmpmx, kcomp)
    AS_ALLOCATE(vi=ent_cod, size=nb_ec)
!
    call dismoi("NOM_MAILLA", chno, "CHAM_NO", repk=mesh)
    call dismoi("NUME_EQUA", chno, "CHAM_NO", repk=numeequa)
    call nueq_chck(numeequa, l_error=.true.)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_point)
!
    if (list_node .ne. ' ') then
        call jelira(list_node, 'LONUTI', nb_point)
        call jeveuo(list_node, 'L', vi=p_list_node)
    else
        call wkvect('&&NUACNO.NOEUD', 'V V I', nb_point, vi=p_list_node)
        do i_pt = 1, nb_point
            p_list_node(i_pt) = i_pt
        end do
    end if
!
    call jelira(chno//'.VALE', 'TYPE', cval=type_scal)
    call jeveuo(chno//'.VALE', 'E', kvale)
    if (type_scal(1:1) .eq. 'R') then
        itype = 1
    else if (type_scal(1:1) .eq. 'C') then
        itype = 2
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(nuage//'.NUAV', 'L', jnuav)
    call jeveuo(nuage//'.NUAI', 'L', jnuai)
    nc = zi(jnuai+2)
!
!     --- SI LE CHAMP EST DECRIT PAR 1 "PRNO" ---
!
    call jeveuo(numeequa//'.NUEQ', 'L', vi=nueq)
    call jenonu(jexnom(numeequa//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(numeequa//'.PRNO', i_ligr_mesh), 'L', iaprno)
    do j = 1, nb_point
        nume_pt = p_list_node(j)
        ival = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+1)
        ncmp = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2)
        if (ncmp .eq. 0) goto 210
        do i_ec = 1, nb_ec
            ent_cod(i_ec) = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2+i_ec)
        end do
        icompt = 0
        do icmp = 1, ncmpmx
            if (exisdg(ent_cod, icmp)) then
                icompt = icompt+1
                ieq = nueq(ival-1+icompt)
                k = nc*(j-1)+icompt
                if (itype .eq. 1) then
                    zr(kvale-1+ieq) = zr(jnuav+k-1)
                else
                    zc(kvale-1+ieq) = zc(jnuav+k-1)
                end if
            end if
        end do
210     continue
    end do
!
    call jedetr('&&NUACNO.NOMCMP')
    AS_DEALLOCATE(vi=ent_cod)
    call jedetr('&&NUACNO.NOEUD')
    call jedema()
end subroutine
