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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine caddlp(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/afddli.h"
#include "asterfort/aflrch.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/char_beam_lcs.h"
#include "asterfort/char_excl_keyw.h"
#include "asterfort/char_read_keyw.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/getnode.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'DDL_POUTRE'
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
    integer(kind=8) :: ddlimp(n_max_keyword)
    real(kind=8) :: valimr(n_max_keyword)
    complex(kind=8) :: valimc(n_max_keyword)
    character(len=8) :: valimf(n_max_keyword)
    character(len=16) :: keywordlist(n_max_keyword)
!
!
    integer(kind=8) :: cmp_nb_glo
    parameter(cmp_nb_glo=6)
    integer(kind=8) :: cmp_acti_glo(cmp_nb_glo)
    real(kind=8) :: cmp_valr_glo(cmp_nb_glo)
    complex(kind=8) :: cmp_valc_glo(cmp_nb_glo)
    character(len=8) :: cmp_valf_glo(cmp_nb_glo)
    character(len=16) :: cmp_name_glo(cmp_nb_glo)
!
    integer(kind=8) :: nddli, iocc, ibid, ino
    integer(kind=8) :: ier, nbec, nbnoeu, n_keyword
    integer(kind=8) :: jdirec, nume_node
    integer(kind=8) :: jnom, nbcmp, jprnm
    real(kind=8) :: zero
    character(len=24) :: keywordexcl
    character(len=4) :: coef_type
    integer(kind=8) :: n_keyexcl
    integer(kind=8) :: i_keyword
    character(len=24) :: list_node
    integer(kind=8) :: jlino
    integer(kind=8) :: nb_node, geomDime
    character(len=8) :: k8bid, nomg, name_node
    character(len=16) :: keywordfact, keyword
    character(len=19) :: lisrel, k19bid, modelLigrel
    character(len=19) :: ncncin
    integer(kind=8), pointer :: dimension(:) => null()
    integer(kind=8), pointer :: icompt(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    keywordfact = 'DDL_POUTRE'
    call getfac(keywordfact, nddli)
    if (nddli .eq. 0) goto 999
!
! - Initializations
!
    lisrel = '&&CADDLP.RLLISTE'
    zero = 0.d0
!
! - Reverse connectivity
!
    ncncin = '&&CADDLP.CONINV'
    call jeexin(ncncin, ier)
    if (ier .eq. 0) call cncinv(mesh, [ibid], 0, 'V', ncncin)
!
    coef_type = 'REEL'
    ASSERT(valeType .eq. 'REEL')

! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call jeveuo(modelLigrel//'.PRNM', 'L', jprnm)
!
! - Create list of excluded keywords for using in load_read_keyw
!
    keywordexcl = '&&CADDLP.KEYWORDEXCL'
    call char_excl_keyw(keywordfact, keywordexcl, n_keyexcl)
!
! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', jnom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmp, k8bid)
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
!
! - Local coordinate system
!
    call jelira(mesh//'.COORDO    .VALE', 'LONMAX', nbnoeu, k8bid)
    nbnoeu = nbnoeu/3
    call wkvect('&&CADDLP.DIRECT', 'V V R', 3*nbnoeu, jdirec)
    AS_ALLOCATE(vi=dimension, size=nbnoeu)
!
! - Loop on factor keyword
!
    do iocc = 1, nddli
!
! ----- Read mesh affectation
!
        list_node = '&&CADDLP.LIST_NODE'
        call getnode(mesh, keywordfact, iocc, 'F', list_node, &
                     nb_node)
        call jeveuo(list_node, 'L', jlino)
!
! ----- Counting components
!
        AS_ALLOCATE(vi=icompt, size=6)
!
! ----- Loop on nodes
!
        do ino = 1, nb_node
            nume_node = zi(jlino-1+ino)
            name_node = int_to_char8(nume_node)
!
! --------- Read affected components and their values
!
            call char_read_keyw(keywordfact, iocc, valeType, n_keyexcl, keywordexcl, &
                                n_max_keyword, n_keyword, keywordlist, ddlimp, valimr, &
                                valimf, valimc)
            ASSERT(n_keyword .le. cmp_nb_glo)
!
! --------- Change components with local coordinate system
!
            call char_beam_lcs(mesh, model, ncncin, keywordfact, iocc, &
                               nume_node, name_node, keywordlist, n_keyword, valimr, cmp_name_glo, &
                               cmp_acti_glo, cmp_valr_glo)
!
! --------- Final linear relation
!
            call afddli(model, geomDime, nbcmp, zk8(jnom), nume_node, name_node, &
                        zi(jprnm-1+(nume_node-1)*nbec+1), dimension(nume_node), &
                        zr(jdirec+3*(nume_node-1)), coef_type, cmp_nb_glo, cmp_name_glo, &
                        cmp_acti_glo, valeType, cmp_valr_glo, cmp_valf_glo, cmp_valc_glo, &
                        icompt, lisrel, .false._1, ibid, ibid, &
                        k19bid, k19bid, k19bid, k19bid, mesh, k19bid)
        end do
        do i_keyword = 1, 6
            keyword = cmp_name_glo(i_keyword)
            if (cmp_acti_glo(i_keyword) .eq. 1) then
                if (icompt(i_keyword) .eq. 0) then
                    call utmess('F', 'CHARGES2_45', sk=keyword)
                end if
            end if
        end do
        AS_DEALLOCATE(vi=icompt)
        call jedetr(list_node)
    end do
!
! - Final linear relation affectation
!
    call aflrch(lisrel, load, 'LIN')
!
    call jedetr('&&CADDLP.DIRECT')
    AS_DEALLOCATE(vi=dimension)
    call jedetr(ncncin)
!
999 continue
!
    call jedema()
!
end subroutine
