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

subroutine xrela_elim(mesh, ds_contact, iden_rela, l_iden_rela, model)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isfonc.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeecra.h"
#include "asterfort/jexnum.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cfdisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
character(len=8), intent(in) :: mesh
character(len=8), intent(in) :: model
type(NL_DS_Contact), intent(in) :: ds_contact
character(len=24), intent(in) :: iden_rela
aster_logical, intent(out) :: l_iden_rela
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! XFEM - Transform linear relations into datastructure for elimination
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  iden_rela        : name of object for identity relations between dof
! In  model            : name of model
! Out l_iden_rela      : .true. if have linear relation to suppress
!
! --------------------------------------------------------------------------------------------------
!
integer :: nbddl1, nbddl2, nbddl3, nbddl4, nb_term_maxi
parameter  (nbddl1=8, nbddl2=36, nbddl3=27, nbddl4=12, nb_term_maxi=12)
character(len=8) :: ddlc1(nbddl1), ddlc2(nbddl2), ddlc3(nbddl3), ddlc4(nbddl4)
integer :: node_nume(nb_term_maxi)
character(len=8) :: node_name(nb_term_maxi), cmp_name(nb_term_maxi)
character(len=8) :: old_node_name, old_cmp_name
character(len=8) :: new_node_name, new_cmp_name, repk
integer :: iret
integer :: nb_crack, nb_dim, nb_edge, nb_iden_rela, nb_iden_term, nb_iden_dof, nb_term
integer :: nb_rela_init, dec, contac, nlag
integer :: i_rela, i_dim, i_edge, i_crack, i_term, i_rela_find, i_rela_idx, i_rela_old
aster_logical :: l_mult_crack, l_rela_find
character(len=14) :: xnrell_crack
character(len=24) :: sdcont_xnrell
character(len=8), pointer :: list_rela(:) => null()
character(len=24), pointer :: v_sdcont_xnrell(:) => null()
integer, pointer :: v_rela_node(:) => null()
integer, pointer :: v_rela_cmp(:) => null()
character(len=8), pointer :: v_sdiden_term(:) => null()
integer, pointer :: v_sdiden_info(:) => null()
integer, pointer :: v_sdiden_dime(:) => null()
integer, pointer :: xfem_cont(:) => null()
!
data ddlc2 /'PRE_FLU','LAG_FLI','LAG_FLS','LAGS_C','LAGS_F1','LAGS_F2',&
            'LAG2_C','LAG2_F1','LAG2_F2','LAG3_C','LAG3_F1','LAG3_F2',&
            'PR2_FLU','LA2_FLI','LA2_FLS','D1X','D1Y','D1Z',&
            'D2X','D2Y','D2Z','D3X','D3Y','D3Z',&
            'PR3_FLU','LA3_FLI','LA3_FLS','V11','V12','V13',&
            'V21','V22','V23','V31','V32','V33'/
data ddlc3 /'PRE_FLU','LAG_FLI','LAG_FLS','LAGS_C','LAGS_F1',&
            'LAG2_C','LAG2_F1','LAG3_C','LAG3_F1',&
            'PR2_FLU','LA2_FLI','LA2_FLS','D1X','D1Y',&
            'D2X','D2Y','D3X','D3Y',&
            'PR3_FLU','LA3_FLI','LA3_FLS','V11','V12',&
            'V21','V22','V31','V32'/
data ddlc4 /'LAGS_C','LAGS_F1','LAGS_F2',&
            'LAG2_C','LAG2_F1','LAG2_F2',&
            'LAG3_C','LAG3_F1','LAG3_F2',&
            'LAG4_C','LAG4_F1','LAG4_F2'/
data ddlc1 /'LAGS_C','LAGS_F1',&
            'LAG2_C','LAG2_F1',&
            'LAG3_C','LAG3_F1',&
            'LAG4_C','LAG4_F1'/
!
! --------------------------------------------------------------------------------------------------
!
nb_rela_init = 0
nb_iden_rela = 0
nb_iden_term = 0
nb_iden_dof  = 0
nb_dim       = cfdisi(ds_contact%sdcont_defi,'NDIM' )
!
! --- TYPE DE CONTACT ET NOMBRE DE MULTIPLICATEURS
!
call jeveuo(model(1:8)//'.XFEM_CONT','L',vi=xfem_cont)
contac = xfem_cont(1)
if(contac.eq.2) nlag = 3
if(contac.eq.1.or.contac.eq.3) nlag = 1
!
! - Get datastructure for linear relations from DEFI_CONTACT
!
    sdcont_xnrell = ds_contact%sdcont_defi(1:16)//'.XNRELL'
    call jeexin(sdcont_xnrell, iret)
    if (iret.eq.0) then
        goto 999
    endif
!
! - HM-XFEM model?
    dec = 0
    call dismoi('EXI_THM', model, 'MODELE', repk=repk)
    if (repk .eq. 'OUI') dec= 3
!
! - Informations about linear relations from DEFI_CONTACT
!
    call jelira(sdcont_xnrell, 'LONMAX', nb_crack)
    call jeveuo(sdcont_xnrell, 'L', vk24 = v_sdcont_xnrell)
!
! - Initial number of linear relations
!
    do i_crack = 1, nb_crack
!
! ----- Current crack
!
        xnrell_crack = v_sdcont_xnrell(i_crack)(1:14)
!
! ----- Number of edges to eliminate
!
        call jeexin(xnrell_crack, iret)
        if (iret.ne.0) then
            call jelira(xnrell_crack, 'LONMAX', nb_edge)
            nb_edge      = nb_edge/2
            nb_rela_init = nb_rela_init + (nb_dim*nlag+dec)*nb_edge
        endif
    end do
!
! - End of treatment if no relations found
!
    if (nb_rela_init .eq. 0) then
        goto 999
    end if
!
! - Create working vector
!
    AS_ALLOCATE(vk8=list_rela, size=nb_rela_init*nb_term_maxi*2)
!
! - Eliminate double

    i_rela = 0
    do i_crack = 1, nb_crack
!
! ----- Current crack
!
        xnrell_crack = v_sdcont_xnrell(i_crack)(1:14)
!
! ----- Number of edges
!
        call jelira(xnrell_crack, 'LONMAX', nb_edge)
        nb_edge = nb_edge/2
!
! ----- Access to nodes
!
        call jeveuo(xnrell_crack, 'L', vi = v_rela_node)
!
! ----- For multi-cracks
!
        call jeexin(xnrell_crack(1:14)//'_LAGR', iret)
        if (iret .eq. 0) then
            l_mult_crack = .false.
        else
            l_mult_crack = .true.
            call jeveuo(xnrell_crack(1:14)//'_LAGR', 'L', vi = v_rela_cmp)
        endif
!
! ----- Loop on edges
!
        do i_edge = 1, nb_edge
!
! --------- Get nodes of linear relation
!
            node_nume(1) = v_rela_node(2*(i_edge-1)+1)
            node_nume(2) = v_rela_node(2*(i_edge-1)+2)
            call jenuno(jexnum(mesh(1:8)//'.NOMNOE', node_nume(1)), node_name(1))
            call jenuno(jexnum(mesh(1:8)//'.NOMNOE', node_nume(2)), node_name(2))
!
! --------- Set linear relation
!
            if (repk.eq.'OUI') then
!
! ------------- Get components name
!
                do i_dim = 1, nlag*nb_dim+3
                    if (l_mult_crack) then
                        if (nb_dim.eq.3) then
                           cmp_name(1) = ddlc2(12*(v_rela_cmp(2*(i_edge-1)+1)-1)+i_dim)
                           cmp_name(2) = ddlc2(12*(v_rela_cmp(2*(i_edge-1)+2)-1)+i_dim)
                        elseif (nb_dim.eq.2) then
                           cmp_name(1) = ddlc3(9*(v_rela_cmp(2*(i_edge-1)+1)-1)+i_dim)
                           cmp_name(2) = ddlc3(9*(v_rela_cmp(2*(i_edge-1)+2)-1)+i_dim)
                        endif
                    else
                        if (nb_dim.eq.3) then
                           cmp_name(1) = ddlc2(i_dim)
                           cmp_name(2) = ddlc2(i_dim)
                        elseif (nb_dim.eq.2) then
                           cmp_name(1) = ddlc3(i_dim)
                           cmp_name(2) = ddlc3(i_dim)
                        endif
                    endif
!
! ------------- Looking for previous relations
!
                    l_rela_find = .false.
                    i_rela_find = 0
                    do i_rela_old = 1, i_rela
                        do i_term = 1, nb_term_maxi
                            old_node_name = list_rela(2*nb_term_maxi*(i_rela_old-1)+2*(i_term-1)+1)
                            old_cmp_name  = list_rela(2*nb_term_maxi*(i_rela_old-1)+2*(i_term-1)+2)
                            if ((old_node_name.eq.node_name(1)).and.&
                                (old_cmp_name.eq.cmp_name(1))) then
                                l_rela_find   = .true.
                                i_rela_find   = i_rela_old
                                new_node_name = node_name(2)
                                new_cmp_name  = cmp_name(2)
                                goto 20
                            endif
                            if ((old_node_name.eq.node_name(2)).and.&
                                (old_cmp_name.eq.cmp_name(2))) then
                                l_rela_find   = .true.
                                i_rela_find   = i_rela_old
                                new_node_name = node_name(1)
                                new_cmp_name  = cmp_name(1)
                                goto 20
                            endif
                        end do
                    end do
!
! ------------- Existing relation
!
 20                 continue
                    if (l_rela_find) then
                        i_rela_idx = 0
                        do i_term = 1, nb_term_maxi
                            node_name(i_term) = &
                                 list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_term-1)+1)
                            if (node_name(i_term).ne.' ') then
                                i_rela_idx = i_rela_idx + 1
                            endif
                        end do
                        if (i_rela_idx.ge.nb_term_maxi) then
                            call utmess('F','XFEM_53')
                        endif
                        i_rela_idx = i_rela_idx + 1
                        list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+1)=new_node_name
                        list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+2)=new_cmp_name
                        nb_iden_term = nb_iden_term + 1
                    endif
!
! ------------- New relation
!
                    if (.not.l_rela_find) then
                        i_rela = i_rela + 1
                        list_rela(2*nb_term_maxi*(i_rela-1)+1) = node_name(1)
                        list_rela(2*nb_term_maxi*(i_rela-1)+2) = cmp_name(1)
                        list_rela(2*nb_term_maxi*(i_rela-1)+3) = node_name(2)
                        list_rela(2*nb_term_maxi*(i_rela-1)+4) = cmp_name(2)
                        nb_iden_term = nb_iden_term + 2
                    endif
                end do
            else
                do i_dim = 1, nlag*nb_dim
!
! ------------- Get components name
!
                    if (l_mult_crack) then
                        cmp_name(1) = ddlc4(3*(v_rela_cmp(2*(i_edge-1)+1)-1)+i_dim)
                        cmp_name(2) = ddlc4(3*(v_rela_cmp(2*(i_edge-1)+2)-1)+i_dim)
                    else
                        if (nb_dim.eq.3) then
                           cmp_name(1) = ddlc4(i_dim)
                           cmp_name(2) = ddlc4(i_dim)
                        elseif (nb_dim.eq.2) then
                           cmp_name(1) = ddlc1(i_dim)
                           cmp_name(2) = ddlc1(i_dim)
                        endif
                    endif
!
! ------------- Looking for previous relations
!
                    l_rela_find = .false.
                    i_rela_find = 0
                    do i_rela_old = 1, i_rela
                        do i_term = 1, nb_term_maxi
                            old_node_name = list_rela(2*nb_term_maxi*(i_rela_old-1)+2*(i_term-1)+1)
                            old_cmp_name  = list_rela(2*nb_term_maxi*(i_rela_old-1)+2*(i_term-1)+2)
                            if ((old_node_name.eq.node_name(1)).and.&
                                (old_cmp_name.eq.cmp_name(1))) then
                                l_rela_find   = .true.
                                i_rela_find   = i_rela_old
                                new_node_name = node_name(2)
                                new_cmp_name  = cmp_name(2)
                                goto 30
                            endif
                            if ((old_node_name.eq.node_name(2)).and.&
                                (old_cmp_name.eq.cmp_name(2))) then
                                l_rela_find   = .true.
                                i_rela_find   = i_rela_old
                                new_node_name = node_name(1)
                                new_cmp_name  = cmp_name(1)
                                goto 30
                            endif
                        end do
                    end do
!
! ------------- Existing relation
!
 30                 continue
                    if (l_rela_find) then
                         i_rela_idx = 0
                         do i_term = 1, nb_term_maxi
                            node_name(i_term) = &
                                 list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_term-1)+1)
                            if (node_name(i_term).ne.' ') then
                                 i_rela_idx = i_rela_idx + 1
                             endif
                        end do
                        if (i_rela_idx.ge.nb_term_maxi) then
                            call utmess('F','XFEM_53')
                        endif
                        old_node_name=list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+1)
                        old_cmp_name =list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+2)
                        if ((old_node_name.ne.new_node_name).and.&
                            (old_cmp_name.ne.new_cmp_name)) then
                            i_rela_idx = i_rela_idx + 1
                            list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+1) = &
                                                                               new_node_name
                            list_rela(2*nb_term_maxi*(i_rela_find-1)+2*(i_rela_idx-1)+2) = &
                                                                                new_cmp_name
                            nb_iden_term = nb_iden_term + 1
                        endif
                    endif
!
! ------------- New relation
!
                    if (.not.l_rela_find) then
                        i_rela = i_rela + 1
                        list_rela(2*nb_term_maxi*(i_rela-1)+1) = node_name(1)
                        list_rela(2*nb_term_maxi*(i_rela-1)+2) = cmp_name(1)
                        list_rela(2*nb_term_maxi*(i_rela-1)+3) = node_name(2)
                        list_rela(2*nb_term_maxi*(i_rela-1)+4) = cmp_name(2)
                        nb_iden_term = nb_iden_term + 2
                    endif
                end do
            endif
        end do
    end do
!
! - Total number of linear relations
!
    nb_iden_rela = i_rela
!
! - Create object for identity relations - Informations
!
    call wkvect(iden_rela(1:19)//'.INFO', 'V V I', 4, vi = v_sdiden_info)
    call wkvect(iden_rela(1:19)//'.DIME', 'V V I', nb_iden_rela, vi = v_sdiden_dime)
!
! - Create object for identity relations - Collection
!
    call jecrec(iden_rela(1:19)//'.COLL', 'V V K8', 'NU', 'CONTIG', 'VARIABLE',&
                nb_iden_rela)
    call jeecra(iden_rela(1:19)//'.COLL', 'LONT', ival=nb_iden_term*2)
!
! - Set objects
!
    do i_rela = 1, nb_iden_rela
!
! ----- Terms of relation
!
        nb_term = 0
        do i_term = 1, nb_term_maxi
            node_name(i_term) = list_rela(2*nb_term_maxi*(i_rela-1)+2*(i_term-1)+1)
            cmp_name(i_term)  = list_rela(2*nb_term_maxi*(i_rela-1)+2*(i_term-1)+2)
            if (node_name(i_term).ne.' ') then
                nb_term = nb_term + 1
            endif
        end do
!
! ----- Create object in collection
!
        call jecroc(jexnum(iden_rela(1:19)//'.COLL', i_rela))
        call jeecra(jexnum(iden_rela(1:19)//'.COLL', i_rela), 'LONMAX',nb_term*2)
        call jeveuo(jexnum(iden_rela(1:19)//'.COLL', i_rela), 'E', vk8 = v_sdiden_term)
!
! ----- Set object in collection
!
        ASSERT(nb_term.ge.2)
        v_sdiden_dime(i_rela) = nb_term
        do i_term = 1, nb_term
            v_sdiden_term(2*(i_term-1)+1) = node_name(i_term)
            v_sdiden_term(2*(i_term-1)+2) = cmp_name(i_term)
        end do
    end do
!
! - Total number of same terms
!
    nb_iden_dof  = nb_iden_term-nb_iden_rela
!
    v_sdiden_info(1) = nb_iden_rela
    v_sdiden_info(2) = nb_iden_term
    v_sdiden_info(3) = nb_iden_dof
!
    AS_DEALLOCATE(vk8=list_rela)
!
999 continue
!
    l_iden_rela = nb_iden_rela.gt.0
!
end subroutine
