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
subroutine nmelco_prep(calc_type, &
                       model, ds_material, ds_contact, &
                       disp_prev, vite_prev, acce_prev, vite_curr, disp_cumu_inst, &
                       nbin, lpain, lchin, &
                       option, time_prev, time_curr, ds_constitutive)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/megeom.h"
!
    character(len=4), intent(in) :: calc_type
    character(len=24), intent(in) :: model
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: disp_prev
    character(len=19), intent(in) :: vite_prev
    character(len=19), intent(in) :: acce_prev
    character(len=19), intent(in) :: vite_curr
    character(len=19), intent(in) :: disp_cumu_inst
    integer(kind=8), intent(in) :: nbin
    character(len=8), intent(out) :: lpain(nbin)
    character(len=19), intent(out) :: lchin(nbin)
    character(len=16), intent(out) :: option
    character(len=19), intent(in) :: time_prev
    character(len=19), intent(in) :: time_curr
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue/LAC methods - Prepare input fields for vector and matrix computation
!
! --------------------------------------------------------------------------------------------------
!
! In  calc_type        : type of computation (vector or matrix)
! In  mesh             : name of mesh
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  ds_contact       : datastructure for contact management
! In  disp_prev        : displacement at beginning of current time
! In  vite_prev        : speed at beginning of current time
! In  vite_curr        : speed at current time
! In  acce_prev        : acceleration at beginning of current time
! In  disp_cumu_inst   : displacement increment from beginning of current time
! In  nbin             : number of input fields
! In  time_prev        : previous time
! In  time_curr        : current time
! In  ds_constitutive  : datastructure for constitutive laws management
! Out lpain            : list of parameters for input fields
! Out lchin            : list of input fields
! Out option           : name of option
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: chgeom, chmlcf, sdappa_psno, sdappa
    aster_logical :: l_cont_cont, l_cont_lac
!
! --------------------------------------------------------------------------------------------------
!
    chgeom = '&&NMELCO.CHGEOM'
    option = ' '
    chmlcf = ' '
    sdappa_psno = ' '
    sdappa = ' '

! - Get contact parameters
    l_cont_cont = cfdisl(ds_contact%sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_lac = cfdisl(ds_contact%sdcont_defi, 'FORMUL_LAC')

! - Select option
    if (calc_type .eq. 'VECT') then
        option = 'CHAR_MECA_CONT'
    elseif (calc_type .eq. 'MATR') then
        option = 'RIGI_CONT'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Geometry field
    call megeom(model, chgeom)

! - <CHELEM> for input field
    if (l_cont_cont) then
        chmlcf = ds_contact%field_input
    end if

! - Special input fields for LAC
    if (l_cont_lac) then
        sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
        chmlcf = ds_contact%field_input
        sdappa_psno = sdappa(1:14)//'.PSNO'
    end if

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PDEPL_M'
    lchin(2) = disp_prev(1:19)
    lpain(3) = 'PDEPL_P'
    lchin(3) = disp_cumu_inst(1:19)
    lpain(4) = 'PVITE_M'
    lchin(4) = vite_prev(1:19)
    lpain(5) = 'PACCE_M'
    lchin(5) = acce_prev(1:19)
    lpain(6) = 'PVITE_P'
    lchin(6) = vite_curr(1:19)
    lpain(7) = 'PMATERC'
    lchin(7) = ds_material%mateco(1:19)
    lpain(8) = 'PCONFR'
    lchin(8) = chmlcf
    lpain(9) = 'PINSTMR'
    lchin(9) = time_prev
    lpain(10) = 'PINSTPR'
    lchin(10) = time_curr
    lpain(11) = 'PCARCRI'
    lchin(11) = ds_constitutive%carcri(1:19)
    lpain(12) = 'PCOMPOR'
    lchin(12) = ds_constitutive%compor(1:19)
    lpain(13) = 'PSNO'
    lchin(13) = sdappa_psno
!
end subroutine
