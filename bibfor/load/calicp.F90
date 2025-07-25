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
subroutine calicp(load, mesh, model, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/indik8.h"
#include "asterfort/aflrch.h"
#include "asterfort/assert.h"
#include "asterfort/char_pair_node.h"
#include "asterfort/dismoi.h"
#include "asterfort/drz12d.h"
#include "asterfort/drz13d.h"
#include "asterfort/getnode.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
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
! Treatment of load LIAISON_COQUE
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
    character(len=16), parameter :: keywordFact = 'LIAISON_COQUE'
    character(len=8) :: k8dummy, nom_noeuds(3)
    character(len=19) :: list_rela, modelLigrel
    integer(kind=8) :: iocc, i_error, icoupl
    integer(kind=8) :: geomDime, nliai, nbec
    character(len=8) :: cmp_name, nomg
    integer(kind=8) :: jnom, nb_cmp
    integer(kind=8) :: cmp_index_dx, cmp_index_dy, cmp_index_dz
    integer(kind=8) :: cmp_index_drx, cmp_index_dry, cmp_index_drz
    character(len=8) :: suffix
    character(len=24) :: list_node_o1, list_node_o2, list_node_i1, list_node_i2
    integer(kind=8) :: nb_node, nb_node_1, nb_node_2
    integer(kind=8) :: j_node_o1, j_node_o2
    integer(kind=8) :: nume_node_1, nume_node_2
    character(len=24) :: list_pair
    integer(kind=8) :: j_list_pair
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordfact, nliai)
    if (nliai .eq. 0) goto 999
!
! - Initializations
!
    list_rela = '&&CALISO.RLLISTE'
    list_node_i1 = '&&CALICP.LIST_NODE_I1'
    list_node_i2 = '&&CALICP.LIST_NODE_I2'
    list_node_o1 = '&&CALICP.LIST_NODE_O1'
    list_node_o2 = '&&CALICP.LIST_NODE_O2'
    list_pair = '&&CALICP.LIST_PAIR'
!
    call wkvect(list_pair, 'V V I', 2, j_list_pair)
!
! - Type
!
    if (valeType .eq. 'COMP') then
        ASSERT(.false.)
    end if
! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
        call utmess('F', 'CHARGES2_6')
    end if

! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', jnom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nb_cmp, k8dummy)
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 10)
!
! - Index in DEPL_R <GRANDEUR> for DX, DY, DZ, DRX, DRY, DRZ
!
    cmp_name = 'DX'
    cmp_index_dx = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    cmp_name = 'DY'
    cmp_index_dy = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    cmp_name = 'DZ'
    cmp_index_dz = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    cmp_name = 'DRX'
    cmp_index_drx = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    cmp_name = 'DRY'
    cmp_index_dry = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    cmp_name = 'DRZ'
    cmp_index_drz = indik8(zk8(jnom), cmp_name, 1, nb_cmp)
    ASSERT(cmp_index_dx .gt. 0)
    ASSERT(cmp_index_dy .gt. 0)
    ASSERT(cmp_index_dz .gt. 0)
    ASSERT(cmp_index_drx .gt. 0)
    ASSERT(cmp_index_dry .gt. 0)
    ASSERT(cmp_index_drz .gt. 0)
!
! - Loop on factor keyword
!
    do iocc = 1, nliai
!
! ----- Read nodes - First list
!
        suffix = '_1'
        call getnode(mesh, keywordfact, iocc, 'F', list_node_i1, &
                     nb_node_1, suffix=suffix)
!
! ----- Read nodes - Second list
!
        suffix = '_2'
        call getnode(mesh, keywordfact, iocc, 'F', list_node_i2, &
                     nb_node_2, suffix=suffix)
!
        if (nb_node_1 .ne. nb_node_2) then
            call utmess('F', 'CHARGES2_8')
        end if
        nb_node = nb_node_1
!
! ----- Create output lists
!
        call wkvect(list_node_o1, 'V V I', nb_node, j_node_o1)
        call wkvect(list_node_o2, 'V V I', nb_node, j_node_o2)
!
! ----- Pairing the two lists with transformation
!
        call char_pair_node(mesh, nb_node, &
                            list_node_i1, list_node_i2, list_node_o1, list_node_o2, i_error)
        if (i_error .ne. 0) then
            call utmess('F', 'CHARGES2_9')
        end if
!
! ----- Compute linear relations
!
        do icoupl = 1, nb_node
            nume_node_1 = zi(j_node_o1-1+icoupl)
            nume_node_2 = zi(j_node_o2-1+icoupl)
            zi(j_list_pair-1+1) = nume_node_1
            zi(j_list_pair-1+2) = nume_node_2
            if (geomDime .eq. 2) then
                call drz12d(mesh, modelLigrel, valeType, 2, list_pair, &
                            cmp_index_drz, list_rela, nom_noeuds)
            else if (geomDime .eq. 3) then
                call drz13d(mesh, modelLigrel, valeType, 2, list_pair, &
                            cmp_index_dx, cmp_index_dy, cmp_index_dz, cmp_index_drx, &
                            cmp_index_dry, cmp_index_drz, list_rela, nom_noeuds)
            end if
        end do
        call jedetr(list_node_i1)
        call jedetr(list_node_i2)
        call jedetr(list_node_o1)
        call jedetr(list_node_o2)
    end do
!
! - Final linear relation affectation
!
    call aflrch(list_rela, load, 'NLIN')
!
    call jedetr(list_pair)
!
999 continue
!
    call jedema()
end subroutine
