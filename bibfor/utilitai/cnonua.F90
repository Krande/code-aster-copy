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

subroutine cnonua(nb_dim, chnoz, list_nodez, nuagez)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/crenua.h"
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
    integer(kind=8), intent(in) :: nb_dim
    character(len=*), intent(in) :: chnoz
    character(len=*), intent(in) :: list_nodez
    character(len=*), intent(in) :: nuagez
!
! --------------------------------------------------------------------------------------------------
!
! Convert CHAM_NO to NUAGE
!
! --------------------------------------------------------------------------------------------------
!
! In  nuage       : name of NUAGE datastructure
! In  chno        : name of CHAM_NO (nodal field) datastructure
! In  nb_dim      : dimension of model
! In  list_node   : list of nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: idx_gd, ncmpmx, nb_ec
    integer(kind=8) :: nb_point, kcoor, kvale, itype
    integer(kind=8) :: nb_cmp_max, i_ec, ianueq, iaprno, nume_pt, ncmp, icompt
    integer(kind=8) :: i_cmp, i_pt, i_dim, i_cmp_mx
    integer(kind=8) :: jnuav, ival, k, ieq, i_ligr_mesh
    character(len=4) :: type_scal
    character(len=8) :: mesh, gran_name
    character(len=19) :: chno, list_node, nuage, numeequa
    aster_logical :: l_crea_nual, prem
    integer(kind=8), pointer :: ent_cod(:) => null()
    integer(kind=8), pointer :: cmp_name(:) => null()
    integer(kind=8), pointer :: p_nuai(:) => null()
    real(kind=8), pointer :: p_nuax(:) => null()
    aster_logical, pointer :: p_nual(:) => null()
    integer(kind=8), pointer :: p_list_node(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    chno = chnoz
    list_node = list_nodez
    nuage = nuagez
    l_crea_nual = .false.
!
    call dismoi("NUM_GD", chno, "CHAM_NO", repi=idx_gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', idx_gd), 'LONMAX', ncmpmx)
    call jenuno(jexnum('&CATA.GD.NOMGD', idx_gd), gran_name)
    nb_ec = nbec(idx_gd)
    AS_ALLOCATE(vi=cmp_name, size=ncmpmx)
    AS_ALLOCATE(vi=ent_cod, size=nb_ec)
!
    call dismoi("NOM_MAILLA", chno, "CHAM_NO", repk=mesh)
    call dismoi("NUME_EQUA", chno, "CHAM_NO", repk=numeequa)
!
    call nueq_chck(numeequa, l_error=.true.)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_point)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', kcoor)
!
    if (list_node .ne. ' ') then
        call jelira(list_node, 'LONUTI', nb_point)
        call jeveuo(list_node, 'L', vi=p_list_node)
    else
        call wkvect('&&CNONUA.NOEUD', 'V V I', nb_point, vi=p_list_node)
        do i_pt = 1, nb_point
            p_list_node(i_pt) = i_pt
        end do
    end if
!
    call jelira(chno//'.VALE', 'TYPE', cval=type_scal)
    call jeveuo(chno//'.VALE', 'L', kvale)
    if (type_scal .eq. 'R') then
        itype = 1
    else if (type_scal .eq. 'C') then
        itype = 2
    else
        ASSERT(.false.)
    end if
!
!
!
!     --- SI LE CHAMP EST DECRIT PAR 1 "PRNO" ---
!     ---------------------------------------------------
    prem = .true.
    call jeveuo(numeequa//'.NUEQ', 'L', ianueq)
    call jenonu(jexnom(numeequa//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(numeequa//'.PRNO', i_ligr_mesh), 'L', iaprno)
    do i_pt = 1, nb_point
        nume_pt = p_list_node(i_pt)
        ncmp = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2)
        if (ncmp .ne. 0) then
            do i_ec = 1, nb_ec
                ent_cod(i_ec) = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2+i_ec)
            end do
            icompt = 0
            do i_cmp = 1, ncmpmx
                if (exisdg(ent_cod, i_cmp)) then
                    icompt = icompt+1
                    cmp_name(i_cmp) = i_cmp
                end if
            end do
            if (prem) then
                nb_cmp_max = icompt
                prem = .false.
            else
                if (nb_cmp_max .ne. icompt) then
                    nb_cmp_max = max(nb_cmp_max, icompt)
                    l_crea_nual = .true.
                end if
            end if
        end if
    end do
!
! - Create NUAGE datastructure
!
    call crenua(nuagez, gran_name, nb_point, nb_dim, nb_cmp_max, &
                l_crea_nual)
!
! - Set .NUAI
!
    call jeveuo(nuage//'.NUAI', 'E', vi=p_nuai)
    p_nuai(1) = nb_point
    p_nuai(2) = nb_dim
    p_nuai(3) = nb_cmp_max
    p_nuai(4) = idx_gd
    p_nuai(5) = itype
    i_cmp = 0
    do i_cmp_mx = 1, ncmpmx
        if (cmp_name(i_cmp_mx) .ne. 0) then
            i_cmp = i_cmp+1
            p_nuai(5+i_cmp) = cmp_name(i_cmp_mx)
        end if
    end do
!
! - Set .NUAX
!
    call jeveuo(nuage//'.NUAX', 'E', vr=p_nuax)
    do i_pt = 1, nb_point
        do i_dim = 1, nb_dim
            p_nuax(nb_dim*(i_pt-1)+i_dim) = zr(kcoor-1+3*(i_pt-1)+i_dim)
        end do
    end do
!
! - Set .NUAV and .NUAL
!
    call jeveuo(nuage//'.NUAV', 'E', jnuav)
    if (l_crea_nual) then
        call jeveuo(nuage//'.NUAL', 'E', vl=p_nual)
    end if
!
!     --- SI LE CHAMP EST DECRIT PAR 1 "PRNO" ---
!
    call jeveuo(numeequa//'.NUEQ', 'L', ianueq)
    call jenonu(jexnom(numeequa//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(numeequa//'.PRNO', i_ligr_mesh), 'L', iaprno)
    do i_pt = 1, nb_point
        nume_pt = p_list_node(i_pt)
        ival = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+1)
        ncmp = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2)
        if (ncmp .ne. 0) then
            do i_ec = 1, nb_ec
                ent_cod(i_ec) = zi(iaprno-1+(nume_pt-1)*(nb_ec+2)+2+i_ec)
            end do
            icompt = 0
            do i_cmp_mx = 1, ncmpmx
                if (exisdg(ent_cod, i_cmp_mx)) then
                    icompt = icompt+1
                    ieq = zi(ianueq-1+ival-1+icompt)
                    k = nb_cmp_max*(i_pt-1)+icompt
                    if (l_crea_nual) then
                        p_nual(k) = .true.
                    end if
                    if (itype .eq. 1) then
                        zr(jnuav+k-1) = zr(kvale-1+ieq)
                    else
                        zc(jnuav+k-1) = zc(kvale-1+ieq)
                    end if
                end if
            end do
        end if
    end do
!
    AS_DEALLOCATE(vi=cmp_name)
    AS_DEALLOCATE(vi=ent_cod)
    call jedetr('&&CNONUA.NOEUD')
    call jedema()
end subroutine
