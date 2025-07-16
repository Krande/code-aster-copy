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
subroutine getelem(mesh, keywordfact, iocc, stop_void, list_elem, &
                   nb_elem, suffix, model, l_keep_propz, l_allz, onAllCells_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: mesh
    character(len=*), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    character(len=1), intent(in) :: stop_void
    integer(kind=8), intent(out) :: nb_elem
    character(len=24), intent(in) :: list_elem
    character(len=8), intent(in), optional :: model
    character(len=*), intent(in), optional :: suffix
    aster_logical, optional, intent(in) :: l_keep_propz
    aster_logical, optional, intent(in) :: l_allz
    aster_logical, optional, intent(out) :: onAllCells_
!
! --------------------------------------------------------------------------------------------------
!
! Read mesh affectation - Elements
!
! --------------------------------------------------------------------------------------------------
!
! Create list of elements:
!  - read MAILLE/GROUP_MA/TOUT keywords
!  - remove by SANS_MAILLE/SANS_GROUP_MA keywords
!  - can use <SUFFIX> to enhance keyword. For instance:
!           suffix = '_1': GROUP_MA -> GROUP_MA_1
!                          MAILLE -> MAILLE_1
!                          SANS_GROUP_MA -> SANS_GROUP_MA_1
!                          SANS_MAILLE -> SANS_MAILLE_1
!           WARNING ->     TOUT -> TOUT
!  - can stop or alarm if no elements in final list
!
! In  mesh         : name of mesh
! In  keywordfact  : factor keyword to read
! In  iocc         : factor keyword index
! In  stop_void    : if nb_elem == 0
!                      'F' - Error
!                      'A' - Error
!                      ' ' - Nothing
! In  list_elem    : list of elements read
! Out nb_elem      : number of elements read
! In  model        : <optional> check elements belongs to model
! In  suffix       : <optional> suffix add to keywords
! IN, OPTIONAL     : L_KEEP_PROP : PRIS EN COMPTE UNIQUEMENT UN MAILLAGE PARALLELE
!    (CELA NE CHANGE RIEN DANS LES AUTRES CAS)
!    POUR UN PARALLEL_MESH, SI TRUE ON NE GARDE QUE LES MAILLES/NOEUDS DONT LE SOUS-DOMAINE
!    EST PROPRIETAIRE SI FALSE ON GARDE TOUT
!    SI L'ARGUMENT N'EST PAS PRESENT ON GARDE TOUT (=FALSE).
! IN, OPTIONAL     : L_ALL : TRUE  : forcer comme TOUT='OUI'
!                            FALSE : par défaut, cas normal
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: moclm(5)
    character(len=16) :: typmcl(5)
    character(len=24) :: list_lect
    integer(kind=8), pointer :: p_list_lect(:) => null()
    character(len=24) :: list_excl
    integer(kind=8), pointer :: p_list_excl(:) => null()
    integer(kind=8), pointer :: p_list_elem(:) => null()
    character(len=24) :: keyword
    character(len=8) :: model_name, suffix_name, answer
    aster_logical :: l_keep_prop, l_all, onAllCells
    integer(kind=8) :: nb_mocl, nb_elem_gl, nocc
    integer(kind=8) :: nb_lect, nb_excl, nb_elim
    integer(kind=8) :: nume_lect, nume_excl
    integer(kind=8) :: i_lect, i_excl, i_elem
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    list_lect = '&&LIST_LECT'
    list_excl = '&&LIST_EXCL'
    nb_elem = 0
    nb_lect = 0
    nb_excl = 0
    model_name = ' '
    suffix_name = ' '
    answer = ' '

    if (present(model)) then
        model_name = model
    end if
    if (present(suffix)) then
        suffix_name = suffix
    end if
    if (present(l_keep_propz)) then
        l_keep_prop = l_keep_propz
    else
        l_keep_prop = ASTER_FALSE
    end if
    if (present(l_allz)) then
        l_all = l_allz
    else
        l_all = ASTER_FALSE
    end if
    onAllCells = ASTER_FALSE
!
! - Read elements
!
    nb_mocl = 0
    nb_mocl = nb_mocl+1
    moclm(nb_mocl) = 'TOUT'
    typmcl(nb_mocl) = 'TOUT'
    keyword = 'GROUP_MA'//suffix_name
    nb_mocl = nb_mocl+1
    moclm(nb_mocl) = keyword
    typmcl(nb_mocl) = 'GROUP_MA'
    keyword = 'MAILLE'//suffix_name
    nb_mocl = nb_mocl+1
    moclm(nb_mocl) = keyword
    typmcl(nb_mocl) = 'MAILLE'
    if (nb_mocl .ne. 0) then
        call reliem(model_name, mesh, 'NU_MAILLE', keywordfact, iocc, &
                    nb_mocl, moclm, typmcl, list_lect, nb_lect, l_keep_prop, l_all)
    end if
    call getvtx(keywordfact, 'TOUT', iocc=iocc, scal=answer, nbret=nocc)
    onAllCells = (answer .eq. "OUI")
!
! - Read elements excludes
!
    nb_mocl = 0
    keyword = 'SANS_GROUP_MA'//suffix_name
    nb_mocl = nb_mocl+1
    moclm(nb_mocl) = keyword
    typmcl(nb_mocl) = 'GROUP_MA'
    keyword = 'SANS_MAILLE'//suffix_name
    nb_mocl = nb_mocl+1
    moclm(nb_mocl) = keyword
    typmcl(nb_mocl) = 'MAILLE'
    if (nb_mocl .ne. 0) then
        call reliem(' ', mesh, 'NU_MAILLE', keywordfact, iocc, &
                    nb_mocl, moclm, typmcl, list_excl, nb_excl)
    end if
!
! - Access to list of elements
!
    if (nb_lect .ne. 0) then
        call jeveuo(list_lect, 'E', vi=p_list_lect)
    end if
!
! - Exclusion of elements in initial list
!
    nb_elim = 0
    if (nb_excl .ne. 0) then
        call jeveuo(list_excl, 'L', vi=p_list_excl)
        do i_excl = 1, nb_excl
            nume_excl = p_list_excl(i_excl)
            do i_lect = 1, nb_lect
                nume_lect = p_list_lect(i_lect)
                if (nume_excl .eq. nume_lect) then
                    nb_elim = nb_elim+1
                    p_list_lect(i_lect) = 0
                end if
            end do
        end do
    end if
    nb_elem = nb_lect-nb_elim
!
! - Final list of elements
!
    i_elem = 0
    if ((nb_elem .ne. 0) .and. (nb_lect .ne. 0)) then
        call wkvect(list_elem, 'V V I', nb_elem, vi=p_list_elem)
        do i_lect = 1, nb_lect
            nume_lect = p_list_lect(i_lect)
            if (nume_lect .ne. 0) then
                i_elem = i_elem+1
                p_list_elem(i_elem) = nume_lect
            end if
        end do
        ASSERT(i_elem .eq. nb_elem)
    end if
!
! - If no elements
!
    nb_elem_gl = nb_elem
    if (isParallelMesh(mesh)) then
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nb_elem_gl)
    end if
!
    if (stop_void .ne. ' ' .and. nb_elem_gl .eq. 0) then
        call utmess(stop_void, 'UTILITY_3', sk=keywordfact)
    end if
!
    if (present(onAllCells_)) then
        onAllCells_ = onAllCells
    end if
!
    call jedetr(list_lect)
    call jedetr(list_excl)
    call jedema()
end subroutine
