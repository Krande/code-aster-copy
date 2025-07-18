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
subroutine lac_crsd(nume_dof, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/vtcreb.h"
!
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! LAC method - Create datastructures for LAC method
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_dof         : name of numbering object (NUME_DDL)
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nt_patch, jv_dummy, nt_elem_slav, vali(2)
    character(len=24) :: sdcont_stat, sdcont_lagc, sdcont_zeta, sdcont_zgpi
    character(len=24) :: sdcont_znmc, sdcont_zpoi, sdcont_zcoe
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'CONTACT5_6')
    end if
!
! - Get parameters
!
    nt_patch = ds_contact%nt_patch
    nt_elem_slav = cfdisi(ds_contact%sdcont_defi, 'NTMAE')
!
! - Print
!
    vali(1) = nt_patch
    vali(2) = nt_elem_slav
    call utmess('I', 'CONTACT5_7', ni=2, vali=vali)
!
! - Create objects
!
    sdcont_stat = ds_contact%sdcont_solv(1:14)//'.STAT'
    sdcont_zeta = ds_contact%sdcont_solv(1:14)//'.ZETA'
    sdcont_zgpi = ds_contact%sdcont_solv(1:14)//'.ZGPI'
    sdcont_zcoe = ds_contact%sdcont_solv(1:14)//'.ZCOE'
    sdcont_znmc = ds_contact%sdcont_solv(1:14)//'.ZNMC'
    sdcont_zpoi = ds_contact%sdcont_solv(1:14)//'.ZPOI'
    sdcont_lagc = ds_contact%sdcont_solv(1:14)//'.LAGC'
    call wkvect(sdcont_stat, 'V V I', nt_patch, jv_dummy)
    call wkvect(sdcont_zeta, 'V V I', nt_patch, jv_dummy)
    call wkvect(sdcont_lagc, 'V V R', nt_patch, jv_dummy)
    call wkvect(sdcont_zgpi, 'V V R', nt_patch, jv_dummy)
    call wkvect(sdcont_zpoi, 'V V R', nt_patch, jv_dummy)
    call wkvect(sdcont_znmc, 'V V I', nt_patch, jv_dummy)
    call wkvect(sdcont_zcoe, 'V V R', nt_patch, jv_dummy)
!
! - Forces to solve
!
    call vtcreb(ds_contact%cneltc, 'V', 'R', nume_ddlz=nume_dof)
    ds_contact%l_cneltc = ASTER_TRUE
!
end subroutine
