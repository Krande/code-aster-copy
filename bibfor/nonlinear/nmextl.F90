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

subroutine nmextl(mesh, model, keyw_fact, i_keyw_fact, field_type, &
                  field_disc, list_node, list_elem, nb_node, nb_elem, &
                  type_extr)
!
    implicit none
!
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/getvtx.h"
#include "asterfort/getelem.h"
#include "asterfort/getnode.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=*), intent(in) :: mesh
    character(len=*), intent(in) :: model
    character(len=16), intent(in) :: keyw_fact
    integer(kind=8), intent(in) :: i_keyw_fact
    character(len=24), intent(in) :: field_type
    character(len=4), intent(in) :: field_disc
    integer(kind=8), intent(out) :: nb_node
    integer(kind=8), intent(out) :: nb_elem
    character(len=24), intent(in) :: list_node
    character(len=24), intent(in) :: list_elem
    character(len=8), intent(out) :: type_extr
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Extraction (OBSERVATION/SUIVI_DDL) utilities
!
! Get topology (nodes or elements) and type of extraction for field
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  keyw_fact        : factor keyword to read extraction parameters
! In  i_keyw_fact      : index of keyword to read extraction parameters
! In  field_type       : type of field (name in results datastructure)
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  list_node        : name of object contains list of nodes
! Out nb_node          : number of nodes
! In  list_elem        : name of object contains list of elements
! Out nb_elem          : number of elements
! Out type_extr        : type of extraction
!                'MIN'      VALEUR MINI SUR TOUTES LES MAILLES/NOEUDS
!                'MAX'      VALEUR MAXI SUR TOUTES LES MAILLES/NOEUDS
!                'MOY'      VALEUR MOYENNE TOUTES LES MAILLES/NOEUDS
!                'MINI_ABS' VALEUR MINI EN ABSOLU SUR TOUTES LES
!                          MAILLES/NOEUDS
!                'MAXI_ABS' VALEUR MAXI EN ABSOLU SUR TOUTES LES
!                          MAILLES/NOEUDS
!                'VALE'     VALEUR TOUTES LES MAILLES/NOEUDS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nocc, nbobj
    aster_logical :: l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    nb_node = 0
    nb_elem = 0
    type_extr = 'VALE'
    l_pmesh = isParallelMesh(mesh)
!
! - Type of extraction on field
!
    call getvtx(keyw_fact, 'EVAL_CHAM', iocc=i_keyw_fact, scal=type_extr, nbret=nocc)
    if (nocc .eq. 0) then
        type_extr = 'VALE'
        call utmess('A', 'EXTRACTION_5', sk=field_type)
    end if
!
! - Get list of nodes
!
    if (field_disc .eq. 'NOEU') then
        call getnode(mesh, keyw_fact, i_keyw_fact, ' ', list_node, &
                     nb_node, model)
        nbobj = nb_node
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbobj)
        end if
        if (nbobj .eq. 0) then
            call utmess('F', 'EXTRACTION_3', sk=field_type)
        end if
    end if
!
! - Get list of elements
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        call getelem(mesh, keyw_fact, i_keyw_fact, ' ', list_elem, &
                     nb_elem, model=model)
        nbobj = nb_elem
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbobj)
        end if
        if (nbobj .eq. 0) then
            call utmess('F', 'EXTRACTION_4', sk=field_type)
        end if
    end if
!
end subroutine
