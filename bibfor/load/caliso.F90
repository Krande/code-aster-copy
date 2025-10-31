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
subroutine caliso(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/indik8.h"
#include "asterfort/aflrch.h"
#include "asterfort/armin.h"
#include "asterfort/assert.h"
#include "asterfort/char_excl_keyw.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/drz12d.h"
#include "asterfort/drz13d.h"
#include "asterfort/exisdg.h"
#include "asterfort/getnode.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/solide_tran.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load LIAISON_SOLIDE
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
    character(len=16), parameter :: keywordFact = 'LIAISON_SOLIDE'
    integer(kind=8) :: iocc
    integer(kind=8) :: jnom, n1
    integer(kind=8) :: i_no
    integer(kind=8) :: nb_cmp, nbec, geomDime, nbocc
    character(len=24) :: list_node
    integer(kind=8) :: jlino, numnoe
    integer(kind=8) :: nb_node
    character(len=8) :: nomg
    real(kind=8) :: dist_mini, dist
    integer(kind=8) :: dim, k
    character(len=1) :: kdim
    character(len=8) :: cmp_name, type_rela, nom_noeuds_tmp(4)
    character(len=8), pointer :: prdso(:) => null()
    integer(kind=8), pointer :: prnm(:) => null()
    integer(kind=8), pointer :: prnm1(:) => null()
    character(len=19) :: list_rela
    character(len=19) :: modelLigrel
    character(len=24) :: keywordexcl
    integer(kind=8) :: n_keyexcl, nuti
    integer(kind=8) :: cmp_index_dx, cmp_index_dy, cmp_index_dz
    integer(kind=8) :: cmp_index_drx, cmp_index_dry, cmp_index_drz
    aster_logical :: l_rota_2d, l_rota_3d
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nbocc)
    if (nbocc .eq. 0) goto 999
!
! - Initializations
!
    list_rela = '&&CALISO.RLLISTE'
!
    l_rota_2d = .false.
    l_rota_3d = .false.

! - Type
    ASSERT(valeType .eq. 'REEL')

! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
        call utmess('F', 'CHARGES2_6')
    end if
    call jeveuo(modelLigrel//'.PRNM', 'L', vi=prnm)

!
! - Minimum distance
!
    dist = armin(mesh)
    ASSERT(dist .gt. 0.d0)
!
! - Create list of excluded keywords for using in char_read_keyw
!
    keywordexcl = '&&CALISO.KEYWORDEXCL'
    call char_excl_keyw(keywordFact, keywordexcl, n_keyexcl)
!
! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', jnom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nb_cmp)
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
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
    do iocc = 1, nbocc
        nom_noeuds_tmp(1:4) = ' '
!
! ----- Minimum distance
!
        call getvr8(keywordFact, 'DIST_MIN', iocc=iocc, scal=dist_mini, nbret=n1)
        if (n1 .eq. 0) dist_mini = dist*1.d-3
!
! ----- Read mesh affectation
!
        list_node = '&&CALISO.LIST_NODE'
        call getnode(mesh, keywordFact, iocc, 'F', list_node, nb_node)
        call jeveuo(list_node, 'L', jlino)
!
! ----- Only one node: nothing to do
!
        if (nb_node .eq. 1) then
            call utmess('I', 'CHARGES2_17')
            cycle
        end if
!
! ----- Model: 2D
!
        if (geomDime .eq. 2) then
!
! --------- Is any node has DRZ dof ?
!
            l_rota_2d = .false.
            do i_no = 1, nb_node
                numnoe = zi(jlino+i_no-1)
                prnm1 => prnm((numnoe-1)*nbec+1:(numnoe-1)*nbec+nbec)
                if (exisdg(prnm1, cmp_index_drz)) then
                    l_rota_2d = .true.
                    goto 40
                end if
            end do
40          continue
!
! --------- Compute linear relations
!
            if (l_rota_2d) then
                call drz12d(mesh, modelLigrel, valeType, nb_node, list_node, &
                            cmp_index_drz, list_rela, nom_noeuds_tmp)
                type_rela = "ROTA2D"
            else
                call solide_tran('2D', mesh, valeType, dist_mini, nb_node, list_node, &
                                 list_rela, nom_noeuds_tmp, dim)
                if (dim .eq. 0) then
                    type_rela = "LIN"
                else
                    call codent(dim, 'D0', kdim)
                    type_rela = "2D"//kdim
                end if
            end if
!
! ----- Model: 3D
!
        else if (geomDime .eq. 3) then
!
! --------- Is any node has rotation dofs ?
!
            l_rota_3d = .false.
            do i_no = 1, nb_node
                numnoe = zi(jlino+i_no-1)
                prnm1 => prnm((numnoe-1)*nbec+1:(numnoe-1)*nbec+nbec)
                if (exisdg(prnm1, cmp_index_drx) .and. &
                    exisdg(prnm1, cmp_index_dry) .and. &
                    exisdg(prnm1, cmp_index_drz)) then
                    l_rota_3d = .true.
                    goto 50
                end if
            end do
50          continue
!
! --------- Compute linear relations
!
            if (l_rota_3d) then
                call drz13d(mesh, modelLigrel, valeType, nb_node, list_node, &
                            cmp_index_dx, cmp_index_dy, cmp_index_dz, cmp_index_drx, &
                            cmp_index_dry, cmp_index_drz, list_rela, nom_noeuds_tmp)
                type_rela = "ROTA3D"
            else
                call solide_tran('3D', mesh, valeType, dist_mini, nb_node, list_node, &
                                 list_rela, nom_noeuds_tmp, dim)
                if (dim .eq. 0) then
                    type_rela = "LIN"
                else
                    call codent(dim, 'D0', kdim)
                    type_rela = "3D"//kdim
                end if
            end if
        end if
!
!       - Final linear relation affectation
!
        call aflrch(list_rela, load, type_rela, elim='NON')
        call detrsd('LISTE_RELA', list_rela)

!       -- remplissage de l'objet .PRDSO :
        call jeveuo(load//'.DUAL.PRDSO', 'E', vk8=prdso)
        call jelira(load//'.DUAL.PRDK', 'LONUTI', ival=nuti)
        do k = 1, 4
            prdso(4*(nuti-1)+k) = nom_noeuds_tmp(k)
        end do
!
        call jedetr(list_node)
!
    end do
!
    call jedetr(keywordexcl)
!
999 continue
    call jedema()
end subroutine
