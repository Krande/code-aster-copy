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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgComporNodes(result, nume_ordre, nb_point, fondNoeudNume, compValues)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calcG_type.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
!
    character(len=8), intent(in) :: result
    integer(kind=8), intent(in) :: nume_ordre, nb_point, fondNoeudNume(*)
    character(len=8), pointer :: compValues(:)
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Create CARTE Comportement for ELAS_INCR
!
! IO compor    : name of <CARTE> COMPORTEMENT
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: connex_inv = '&&CGCMPS.CONINV'
    character(len=8) :: model, mesh
    character(len=16) :: rela_comp
    character(len=19) :: ligrmo, compor
    integer(kind=8) :: nb_vale, nb_zone, nb_cmp_max, iret
    integer(kind=8) :: i_node, i_cell, nb_node_cell, node_nume, cell_nume
    integer(kind=8) :: compor_node, compor_cell, zone_nume
    integer(kind=8), pointer :: comporPtma(:) => null()
    integer(kind=8), pointer :: v_coninv(:) => null()
    integer(kind=8), pointer :: v_coninv_longcum(:) => null()
    integer(kind=8), pointer :: v_compor_desc(:) => null()
    character(len=16), pointer :: v_compor_vale(:) => null()
!
    call jemarq()
!
    call dismoi('MODELE', result, 'RESULTAT', repk=model)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    compValues(:) = '        '
!
! --- Read COMPOR <CARTE> in RESULT
!
    call rsexch(' ', result, 'COMPORTEMENT', nume_ordre, compor, iret)
!
    if (iret .ne. 0) then
! - Everyting is elastic
        compValues(1:nb_point) = "ELAS"
        go to 999
    end if
!
! - Prepare COMPOR field
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
    call etenca(compor, ligrmo, iret)
    ASSERT(iret == 0)
    call jeveuo(compor//'.PTMA', 'L', vi=comporPtma)
!
! - Access to COMPOR
!
    call jeveuo(compor//'.DESC', 'L', vi=v_compor_desc)
    call jeveuo(compor//'.VALE', 'E', vk16=v_compor_vale)
    call jelira(compor//'.VALE', 'LONMAX', nb_vale)
    nb_zone = v_compor_desc(3)
    nb_cmp_max = nb_vale/v_compor_desc(2)
!
! ----- Inverse connectivity and mesh paramters
!
    call cncinv(mesh, [0], 0, 'V', connex_inv)
!
    do i_node = 1, nb_point
        compor_node = -1
        node_nume = fondNoeudNume(i_node)
! --- For nb_point_fond: behavior not identified (could be possible but harder)
        if (node_nume <= 0) then
            compValues(i_node) = "XXXXXXXX"
            cycle
        end if
! --------- Get elements attached to current node
        call jeveuo(jexatr(connex_inv, 'LONCUM'), 'L', vi=v_coninv_longcum)
        nb_node_cell = v_coninv_longcum(node_nume+1)-v_coninv_longcum(node_nume)
        call jeveuo(jexnum(connex_inv, node_nume), 'L', vi=v_coninv)
! --------- Loop on elements attached to current node
        do i_cell = 1, nb_node_cell
            cell_nume = v_coninv(i_cell)
            zone_nume = comporPtma(cell_nume)
            rela_comp = v_compor_vale(nb_cmp_max*(zone_nume-1)+RELA_NAME)
            if (rela_comp == "ELAS" .or. rela_comp == "ELAS_FO") then
                compor_cell = ELAS
            elseif (rela_comp(1:10) == "ELAS_VMIS_") then
                compor_cell = ELAS_NL
            elseif (rela_comp(1:10) == "VMIS_ISOT_") then
                compor_cell = PLAS
            else
                ASSERT(ASTER_FALSE)
            end if
            compor_node = max(compor_node, compor_cell)
        end do
!
        select case (compor_node)
        case (ELAS)
            compValues(i_node) = "ELAS"
        case (ELAS_NL)
            compValues(i_node) = "ELAS_NL"
        case (PLAS)
            compValues(i_node) = "PLAS"
        case default
            ASSERT(ASTER_FALSE)
        end select
    end do
!
999 continue
!
    call jedema()
!
end subroutine
