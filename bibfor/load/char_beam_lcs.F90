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

subroutine char_beam_lcs(mesh, model, connex_inv, keywordfact, iocc, &
                         node_nume, node_name, cmp_name_loc, n_keyword, cmp_valr_loc, &
                         cmp_name_glo, cmp_acti_glo, cmp_valr_glo)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/matloc.h"
#include "asterfort/reliem.h"
#include "asterfort/utpvlg.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: connex_inv
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    integer(kind=8), intent(in) :: node_nume
    integer(kind=8), intent(in) :: n_keyword
    character(len=8), intent(in) :: node_name
    character(len=16), intent(in) :: cmp_name_loc(6)
    real(kind=8), intent(in) :: cmp_valr_loc(6)
    character(len=16), intent(out) :: cmp_name_glo(6)
    integer(kind=8), intent(out) :: cmp_acti_glo(6)
    real(kind=8), intent(out) :: cmp_valr_glo(6)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Change components with local coordinate system for beams at node
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh           : meshing
! In  model          : model
! In  connex_inv     : inverse connectivity of mesh (nodes -> elements)
! In  keywordfact    : factor cmp_name to read
! In  iocc           : factor cmp_name index in AFFE_CHAR_MECA
! In  node_name      : name of node
! In  node_nume      : number of node
! In  cmp_name_loc   : list of components in local coordinate system
! In  cmp_valr_loc   : values (if real) of components in local coordinate system
! Out cmp_name_glo   : list of components in global coordinate system
! Out cmp_acti_glo   : 1 if component affected, 0 else
! Out cmp_valr_glo   : values (if real) of components in global coordinate system
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: matr_glob_loca(3, 3)
    real(kind=8) :: rln1(3), rgn1(3)
    real(kind=8) :: dloc(3), dglo(3)
    integer(kind=8) :: i_direc, i_cmp
    character(len=16) :: keyw_name(2), keyw_type(2)
    character(len=16) :: cmp_name, list_cmp(6)
    character(len=24) :: list_repe_elem
    integer(kind=8) :: nb_repe_elem, j_repe_elem
!
    data list_cmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    ASSERT(n_keyword .le. 6)
!
! - Initializations
!
    do i_cmp = 1, 6
        cmp_valr_glo(i_cmp) = 0.d0
        cmp_acti_glo(i_cmp) = 0
        cmp_name_glo(i_cmp) = list_cmp(i_cmp)
    end do
!
! - Mesh for local coordinate system
!
    keyw_name(1) = 'MAILLE_REPE'
    keyw_type(1) = 'MAILLE'
    keyw_name(2) = 'GROUP_MA_REPE'
    keyw_type(2) = 'GROUP_MA'
    list_repe_elem = '&&REPE.MAILLE'
    call reliem(model, mesh, 'NU_MAILLE', keywordfact, iocc, &
                2, keyw_name, keyw_type, list_repe_elem, nb_repe_elem)
    if (nb_repe_elem .ne. 0) then
        call jeveuo(list_repe_elem, 'L', j_repe_elem)
    else
        j_repe_elem = 1
    end if
!
! - Local coordinate system
!
    call matloc(mesh, connex_inv, keywordfact, iocc, node_nume, &
                node_name, nb_repe_elem, zi(j_repe_elem), matr_glob_loca)
!
! - Translation
!
    do i_direc = 1, 3
        dloc(i_direc) = 0.d0
        rln1(i_direc) = 0.d0
    end do
    do i_cmp = 1, n_keyword
        cmp_name = cmp_name_loc(i_cmp)
        if (cmp_name .eq. 'DX') then
            rln1(1) = 1.d0
            dloc(1) = cmp_valr_loc(i_cmp)
        end if
        if (cmp_name .eq. 'DY') then
            rln1(2) = 1.d0
            dloc(2) = cmp_valr_loc(i_cmp)
        end if
        if (cmp_name .eq. 'DZ') then
            rln1(3) = 1.d0
            dloc(3) = cmp_valr_loc(i_cmp)
        end if
    end do
!
    call utpvlg(1, 3, matr_glob_loca, dloc, dglo)
    call utpvlg(1, 3, matr_glob_loca, rln1, rgn1)
!
    if (rgn1(1) .ne. 0.d0) then
        cmp_valr_glo(1) = dglo(1)
        cmp_acti_glo(1) = 1
    end if
    if (rgn1(2) .ne. 0.d0) then
        cmp_valr_glo(2) = dglo(2)
        cmp_acti_glo(2) = 1
    end if
    if (rgn1(3) .ne. 0.d0) then
        cmp_valr_glo(3) = dglo(3)
        cmp_acti_glo(3) = 1
    end if
!
! - Rotation
!
    do i_direc = 1, 3
        dloc(i_direc) = 0.d0
        rln1(i_direc) = 0.d0
    end do
    do i_cmp = 1, n_keyword
        cmp_name = cmp_name_loc(i_cmp)
        if (cmp_name .eq. 'DRX') then
            rln1(1) = 1.d0
            dloc(1) = cmp_valr_loc(i_cmp)
        end if
        if (cmp_name .eq. 'DRY') then
            rln1(2) = 1.d0
            dloc(2) = cmp_valr_loc(i_cmp)
        end if
        if (cmp_name .eq. 'DRZ') then
            rln1(3) = 1.d0
            dloc(3) = cmp_valr_loc(i_cmp)
        end if
    end do
!
    call utpvlg(1, 3, matr_glob_loca, dloc, dglo)
    call utpvlg(1, 3, matr_glob_loca, rln1, rgn1)
!
    if (rgn1(1) .ne. 0.d0) then
        cmp_valr_glo(4) = dglo(1)
        cmp_acti_glo(4) = 1
    end if
    if (rgn1(2) .ne. 0.d0) then
        cmp_valr_glo(5) = dglo(2)
        cmp_acti_glo(5) = 1
    end if
    if (rgn1(3) .ne. 0.d0) then
        cmp_valr_glo(6) = dglo(3)
        cmp_acti_glo(6) = 1
    end if
!
    call jedema()
end subroutine
