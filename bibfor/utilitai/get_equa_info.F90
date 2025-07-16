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

subroutine get_equa_info(nume_ddlz, i_equa, type_equa, nume_nodez, nume_cmpz, &
                         nume_cmp_lagrz, nume_subsz, nume_linkz, nb_node_lagr, list_node_lagr, &
                         ligrelz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeexin.h"
#include "asterfort/get_lagr_info.h"
!
!
    character(len=*), intent(in) :: nume_ddlz
    integer(kind=8), intent(in) :: i_equa
    character(len=*), intent(out) :: type_equa
    integer(kind=8), optional, intent(out) :: nume_nodez
    integer(kind=8), optional, intent(out) :: nume_cmpz
    integer(kind=8), optional, intent(out) :: nume_subsz
    integer(kind=8), optional, intent(out) :: nume_linkz
    integer(kind=8), optional, intent(out) :: nume_cmp_lagrz
    integer(kind=8), optional, intent(out) :: nb_node_lagr
    integer(kind=8), optional, pointer :: list_node_lagr(:)
    character(len=*), optional, intent(out) :: ligrelz
!
! --------------------------------------------------------------------------------------------------
!
! Get information about dof (node, component, etc.)
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_ddl       : name of numbering (NUME_DDL)
! In  i_equa         : index of equation
! Out type_equa      : type of dof
!                 / 'A' : physical dof (node+component)
!                 / 'B' : Lagrange dof (boundary condition) simple given boundary condition
!                 / 'C' : Lagrange dof (boundary condition) linear relation
!                 / 'D' : generalized dof - Substructuring
!                 / 'E' : generalized dof - Links
! Out nume_node      : global node index in mesh
! Out nume_cmp       : global component index in GRANDEUR
! Out nume_cmp_lagr  : global component index in GRANDEUR for node linked to Lagrange node
! Out nume_subs      : index of substructure (generalized dof)
! Out nume_link      : index of kinematic link (generalized dof)
! Out nb_node_lagr   : number of nodes linked to lagrange dof
! Out list_node_lagr : pointer to list of nodes linked to lagrange dof
! Out ligrel         : name of LIGREL for non-physical node (Lagrange)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: lili, nueq, orig, desc, deeq
    integer(kind=8) :: nume_node, nume_cmp, nume_subs, nume_link, nume_cmp_lagr
    character(len=19) :: nume_equa, ligrel
    character(len=14) :: nume_ddl
    integer(kind=8) :: iexi, isst
    integer(kind=8) :: idx_gd
    integer(kind=8) :: ino, icmp
    logical :: l_gene
    integer(kind=8), pointer :: p_nueq(:) => null()
    integer(kind=8), pointer :: p_desc(:) => null()
    integer(kind=8), pointer :: p_deeq(:) => null()
    integer(kind=8), pointer :: p_orig(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_ddl = nume_ddlz
    type_equa = '?'
    nume_node = 0
    nume_cmp = 0
    nume_subs = 0
    nume_link = 0
    nume_cmp_lagr = 0
!
! - Get name of nume_equa
!
    call dismoi('NUME_EQUA', nume_ddl, 'NUME_DDL', repk=nume_equa)
!
! - NUME_EQUA or NUME_EQUA_GENE ?
!
    call jeexin(nume_equa//'.DESC', iexi)
    l_gene = (iexi .gt. 0)
!
! - Objects in NUME_EQUA/NUME_EQUA_GENE
!
    deeq = nume_equa(1:19)//'.DEEQ'
    lili = nume_equa(1:19)//'.LILI'
    desc = nume_equa(1:19)//'.DESC'
    orig = nume_equa(1:19)//'.ORIG'
    nueq = nume_equa(1:19)//'.NUEQ'
!
    call jeveuo(deeq, 'L', vi=p_deeq)
    call jeveuo(nueq, 'L', vi=p_nueq)
!
    if (l_gene) then
        call jeveuo(desc, 'L', vi=p_desc)
        ASSERT(p_desc(1) .eq. 2)
        isst = p_deeq(2*(i_equa-1)+2)
        if (isst .gt. 0) then
            type_equa = 'D'
            call jeexin(jexnum(orig, 1), iexi)
            if (iexi .gt. 0) then
                call jeveuo(jexnum(orig, 1), 'L', vi=p_orig)
                nume_subs = p_orig(isst)
            else
                nume_subs = 0
            end if
        else
            type_equa = 'E'
            isst = -isst
            call jeveuo(jexnum(orig, 2), 'L', vi=p_orig)
            nume_link = p_orig(isst)
        end if
    else
        call dismoi('NUM_GD_SI', nume_ddl, 'NUME_DDL', repi=idx_gd)
!
        ino = p_deeq(2*(i_equa-1)+1)
        icmp = p_deeq(2*(i_equa-1)+2)
!
! ----- Physical node
!
        if (ino .gt. 0 .and. icmp .gt. 0) then
            type_equa = 'A'
            nume_node = ino
            nume_cmp = icmp
            goto 70
        end if
!
! ----- Non-Physical node (Lagrange)
!
        if (ino .gt. 0 .and. icmp .lt. 0) then
            type_equa = 'B'
            call get_lagr_info(nume_equa, i_equa, idx_gd, nb_node_lagr, list_node_lagr, &
                               nume_cmp)
            ASSERT(nb_node_lagr .eq. 1)
            nume_node = list_node_lagr(1)
            nume_cmp_lagr = nume_cmp
            goto 70
        end if
!
! ----- Non-Physical node (Lagrange) - LIAISON_DDL
!
        if (ino .eq. 0 .and. icmp .eq. 0) then
            type_equa = 'C'
            call get_lagr_info(nume_equa, i_equa, idx_gd, nb_node_lagr, list_node_lagr, &
                               ligrelz=ligrel)
            goto 70
        end if
    end if
!
70  continue
!
    if (present(nume_nodez)) then
        nume_nodez = nume_node
    end if
    if (present(nume_cmpz)) then
        nume_cmpz = nume_cmp
    end if
    if (present(nume_cmp_lagrz)) then
        nume_cmp_lagrz = nume_cmp_lagr
    end if
    if (present(nume_subsz)) then
        nume_subsz = nume_subs
    end if
    if (present(nume_linkz)) then
        nume_linkz = nume_link
    end if
    if (present(ligrelz)) then
        ligrelz = ligrel
    end if
!
end subroutine
