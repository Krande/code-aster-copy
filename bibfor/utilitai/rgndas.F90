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
subroutine rgndas(nume_ddlz, i_equa, l_print, type_equaz, name_nodez, &
                  name_cmpz, ligrelz)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/get_equa_info.h"
#include "asterfort/equa_print.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*), intent(in) :: nume_ddlz
    integer(kind=8), intent(in) :: i_equa
    logical, intent(in) :: l_print
    character(len=1), optional, intent(out) :: type_equaz
    character(len=*), optional, intent(out) :: name_nodez
    character(len=*), optional, intent(out) :: name_cmpz
    character(len=*), optional, intent(out) :: ligrelz
!
! --------------------------------------------------------------------------------------------------
!
! Get and/or print information about dof (node, component, etc.)
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_ddl      : name of numbering (NUME_DDL)
! In  i_equa        : index of equation
! In  l_print       : .true. to print equation information
! Out type_equa      : type of dof
!                 / 'A' : physical dof (node+component)
!                 / 'B' : Lagrange dof (boundary condition) simple given boundary condition
!                 / 'C' : Lagrange dof (boundary condition) linear relation
!                 / 'D' : generalized dof - Substructuring
!                 / 'E' : generalized dof - Links
! Out name_node      : name of the node
! Out name_cmp       : name of the component
! Out ligrel         : name of LIGREL for non-physical node (Lagrange)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: idx_gd, iexi
    character(len=1) :: type_equa
    character(len=8) :: mesh, modl_gene, ligrel
    integer(kind=8) :: nume_node, nume_cmp, nume_cmp_lagr, nume_subs, nume_link
    character(len=19) :: nume_equa_gene
    character(len=14) :: nume_ddl
    character(len=8) :: name_node, name_cmp, name_cmp_lagr, name_subs
    character(len=8), pointer :: p_cata_nomcmp(:) => null()
    character(len=24), pointer :: p_refe(:) => null()
    integer(kind=8) :: nb_node_lagr
    integer(kind=8), pointer:: list_node_lagr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_ddl = nume_ddlz
    nume_equa_gene = nume_ddl(1:14)//'.NUME'
    ligrel = ' '
    name_node = ' '
    name_cmp = ' '
    name_cmp_lagr = ' '
    name_subs = ' '
    nb_node_lagr = 0

    call dismoi('NOM_MAILLA', nume_ddlz, 'NUME_DDL', repk=mesh)
    call dismoi('NUM_GD_SI', nume_ddlz, 'NUME_DDL', repi=idx_gd)
!
! - No information about mesh (discrete contact case)
!
    if (mesh .eq. ' ' .or. idx_gd .eq. 0) then
        type_equa = 'A'
        name_node = 'Unknown'
        name_cmp = 'Unknown'
        goto 99
    end if
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', idx_gd), 'L', vk8=p_cata_nomcmp)
!
! - Get information about dof
!
    call get_equa_info(nume_ddlz, i_equa, type_equa, nume_node, nume_cmp, &
                       nume_cmp_lagr, nume_subs, nume_link, nb_node_lagr, list_node_lagr, &
                       ligrel)
!
! - Physical dof
!
    if (type_equa .eq. 'A') then
        name_node = int_to_char8(nume_node)
        name_cmp = p_cata_nomcmp(nume_cmp)
    end if
!
! - Non-Physical dof (Lagrange)
!
    if (type_equa .eq. 'B') then
        name_node = int_to_char8(nume_node)
        name_cmp = p_cata_nomcmp(nume_cmp)
        ASSERT(name_cmp .eq. 'LAGR')
        name_cmp_lagr = p_cata_nomcmp(nume_cmp_lagr)
    end if
!
! - Non-Physical dof (Lagrange) - LIAISON_DDL
!
    if (type_equa .eq. 'C') then
!
    end if
!
! - Generalized dof - Substructuring
!
    if (type_equa .eq. 'D') then
        name_cmp = 'GEN'
        call jeexin(nume_equa_gene//'.REFE', iexi)
        if (iexi .gt. 0) then
            call jeveuo(nume_equa_gene//'.REFE', 'L', vk24=p_refe)
        else
            call jeveuo(nume_equa_gene//'.REFN', 'L', vk24=p_refe)
        end if
        modl_gene = p_refe(1) (1:8)
        call jeexin(modl_gene//'      .MODG.SSNO', iexi)
        if (iexi .gt. 0) then
            if (nume_subs .eq. 0) then
                name_subs = 'N/A'
            else
                call jenuno(jexnum(modl_gene//'      .MODG.SSNO', nume_subs), name_subs)
            end if
        else
            name_subs = 'N/A'
        end if
        name_node = name_subs
    end if
!
! - Generalized dof - Kinematic link
!
    if (type_equa .eq. 'E') then
        name_node = 'TAR'
        name_cmp = 'LAG'
    end if
!
! - Print equation
!
    if (l_print) then
        call equa_print(mesh, i_equa, type_equa, name_node, name_cmp, &
                        name_cmp_lagr, name_subs, nume_link, nb_node_lagr, list_node_lagr, &
                        ligrel)
    end if
!
99  continue
    if (present(name_nodez)) then
        name_nodez = name_node
    end if
    if (present(name_cmpz)) then
        name_cmpz = name_cmp
    end if
    if (present(type_equaz)) then
        type_equaz = type_equa
    end if
    if (present(ligrelz)) then
        ligrelz = ligrel
    end if
!
    if (nb_node_lagr .gt. 0) then
        AS_DEALLOCATE(vi=list_node_lagr)
    end if
!
end subroutine
