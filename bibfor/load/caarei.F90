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
!
subroutine caarei(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/afddli.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/catang.h"
#include "asterfort/char_excl_keyw.h"
#include "asterfort/char_read_keyw.h"
#include "asterfort/char_read_val.h"
#include "asterfort/char_xfem.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/getnode.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xddlim.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'ARETE_IMPO'
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh        : mesh
! In  load        : load
! In  model       : model
! In  valeType    : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: n_max_keyword = 320
    integer(kind=8) :: nbterm(n_max_keyword)
    real(kind=8) :: vale_real(n_max_keyword)
    complex(kind=8) :: vale_cplx(n_max_keyword)
    character(len=8) :: vale_func(n_max_keyword)
    character(len=16) :: keywordlist(n_max_keyword)
!
    integer(kind=8) :: jdirec, jprnm
    integer(kind=8) :: iocc, nume_node, ibid
    integer(kind=8) :: ino, inom, idim
    integer(kind=8) :: nbnoeu, narei, nbcmp, nbec, geomDime
    real(kind=8) :: repe_defi(3)
    integer(kind=8) :: repe_type
    real(kind=8) :: coef_real_unit
    complex(kind=8) :: coef_cplx_unit
    integer(kind=8) :: i_keyword, n_keyword
    character(len=24) :: list_node, list_elem
    integer(kind=8) :: jlino, jlima
    integer(kind=8) :: nb_node, nb_elem
    character(len=4) :: coef_type
    character(len=8) :: nomg
    character(len=8) :: name_node, dof_name, k8bid
    character(len=16) :: keywordfact, keyword
    character(len=19) :: list_rela
    character(len=19) :: connex_inv, modelLigrel
    character(len=19) :: ch_xfem_stat, ch_xfem_node, ch_xfem_lnno, ch_xfem_ltno, ch_xfem_heav
    integer(kind=8) :: jnoxfl, jnoxfv
    aster_logical :: lxfem, l_ocmp
    aster_logical :: l_dtan
    integer(kind=8) :: val_nb_dtan
    real(kind=8) :: val_r_dtan
    character(len=8) :: val_f_dtan
    complex(kind=8) :: val_c_dtan
    character(len=16) :: val_t_dtan
    character(len=24) :: keywordexcl
    integer(kind=8) :: n_keyexcl
    integer(kind=8), pointer :: icompt(:) => null()
    real(kind=8), pointer :: tangent(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    keywordfact = 'ARETE_IMPO'
    call getfac(keywordfact, narei)
    if (narei .eq. 0) goto 999
!
! - Initializations
!
    list_rela = '&&CAAREI.RLLISTE'
    coef_cplx_unit = (1.d0, 0.d0)
    coef_real_unit = 1.d0
    dof_name = 'DEPL'

! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    if (geomDime .ne. 3) then
        call utmess('F', 'CHARGES2_7', si=geomDime)
    end if
    call jeveuo(modelLigrel//'.PRNM', 'L', jprnm)
!
! - Type of coefficients
!
    coef_type = 'REEL'
    ASSERT(valeType .eq. 'REEL')
!
! - Create list of excluded keywords
!
    keywordexcl = '&&CAAREI.KEYWORDEXCL'
    call char_excl_keyw(keywordfact, keywordexcl, n_keyexcl)
!
! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', inom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmp)
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
!
! - Local coordinate system (dummy)
!
    call jelira(mesh//'.COORDO    .VALE', 'LONMAX', nbnoeu, k8bid)
    nbnoeu = nbnoeu/3
    call wkvect('&&CAAREI.REPE_DEFI', 'V V R', 3*nbnoeu, jdirec)
!
! - Xfem fields
!
    call char_xfem(mesh, model, lxfem, connex_inv, ch_xfem_stat, &
                   ch_xfem_node, ch_xfem_lnno, ch_xfem_ltno, ch_xfem_heav)
    if (lxfem) then
        call jeveuo(ch_xfem_node//'.CNSL', 'L', jnoxfl)
        call jeveuo(ch_xfem_node//'.CNSV', 'L', jnoxfv)
    end if
!
! - Loop on factor keyword
!
    do iocc = 1, narei
!
! ----- Read mesh affectation
!
        list_node = '&&CAAREI.LIST_NODE'
        list_elem = '&&CAAREI.LIST_ELEM'
        call getnode(mesh, keywordfact, iocc, 'F', list_node, &
                     nb_node)
        call getelem(mesh, keywordfact, iocc, 'F', list_elem, &
                     nb_elem)
        call jeveuo(list_node, 'L', jlino)
        call jeveuo(list_elem, 'L', jlima)
!
! ----- Read keywords and their values except for affectation
!
        call char_read_keyw(keywordfact, iocc, valeType, n_keyexcl, keywordexcl, &
                            n_max_keyword, n_keyword, keywordlist, nbterm, vale_real, &
                            vale_func, vale_cplx)
!
! ----- Detection of DTAN and others
!
        call char_read_val(keywordfact, iocc, 'DTAN', valeType, val_nb_dtan, &
                           val_r_dtan, val_f_dtan, val_c_dtan, val_t_dtan)
        l_dtan = val_nb_dtan .gt. 0
        l_ocmp = n_keyword .gt. 0
!
! ----- Tangents
!
        if (l_dtan) then
            call catang(mesh, nb_elem, zi(jlima), nb_node, zi(jlino))
            call jeveuo('&&CATANG.TANGENT', 'L', vr=tangent)
        end if
!
! ----- If DTAN exists
!
        if (l_dtan) then
            do ino = 1, nb_node
                nume_node = zi(jlino+ino-1)
                name_node = int_to_char8(nume_node)
                do idim = 1, geomDime
                    repe_defi(idim) = tangent(geomDime*(ino-1)+idim)
                end do
!
                if (lxfem) then
                    if (zl(jnoxfl-1+2*nume_node)) then
                        call xddlim(model, dof_name, name_node, nume_node, val_r_dtan, &
                                    val_c_dtan, val_f_dtan, valeType, ibid, list_rela, &
                                    geomDime, repe_defi, jnoxfv, ch_xfem_stat, ch_xfem_lnno, &
                                    ch_xfem_ltno, connex_inv, ch_xfem_heav)
                        goto 115
                    end if
                end if
!
                repe_type = geomDime
                call afrela([coef_real_unit], [coef_cplx_unit], dof_name, name_node, [repe_type], &
                            repe_defi, val_nb_dtan, val_r_dtan, val_c_dtan, val_f_dtan, &
                            coef_type, valeType, 0.d0, list_rela)
!
115             continue
            end do
        end if
!
! ----- If other components exist
!
        if (l_ocmp) then
!
! --------- Counting components
!
            AS_ALLOCATE(vi=icompt, size=n_keyword)
!
! --------- Linear relation
!
            do ino = 1, nb_node
                nume_node = zi(jlino-1+ino)
                name_node = int_to_char8(nume_node)
                call afddli(model, geomDime, nbcmp, zk8(inom), nume_node, name_node, &
                            zi(jprnm-1+(nume_node-1)*nbec+1), 0, zr(jdirec+3*(nume_node-1)), &
                            coef_type, n_keyword, keywordlist, nbterm, valeType, &
                            vale_real, vale_func, vale_cplx, icompt, list_rela, &
                            lxfem, jnoxfl, jnoxfv, ch_xfem_stat, ch_xfem_lnno, &
                            ch_xfem_ltno, connex_inv, mesh, ch_xfem_heav)
!
            end do
!
! --------- Components doesn't exist on all nodes
!
            do i_keyword = 1, n_keyword
                keyword = keywordlist(i_keyword)
                if (icompt(i_keyword) .eq. 0) then
                    call utmess('F', 'CHARGES2_45', sk=keyword)
                end if
            end do
!
            AS_DEALLOCATE(vi=icompt)
!
        end if
!
        call jedetr(list_node)
        call jedetr(list_elem)
!
    end do
!
! - Final linear relation affectation
!
    call aflrch(list_rela, load, 'LIN')
!
    call jedetr('&&CATANG.TANGENT')
    call jedetr('&&CAAREI.REPE_DEFI')
    call jedetr(keywordexcl)
    if (lxfem) then
        call jedetr(connex_inv)
        call detrsd('CHAM_NO_S', ch_xfem_node)
        call detrsd('CHAM_ELEM_S', ch_xfem_stat)
        call detrsd('CHAM_ELEM_S', ch_xfem_lnno)
        call detrsd('CHAM_ELEM_S', ch_xfem_ltno)
    end if
!
999 continue
    call jedema()
end subroutine
