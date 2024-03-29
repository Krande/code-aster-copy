! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine xminit(mesh, model, ds_contact, nume_inst, ds_measure, &
                  sddyna, hat_valinc, list_func_acti)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/copisd.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/xmelem.h"
#include "asterfort/mmbouc.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    type(NL_DS_Contact), intent(inout) :: ds_contact
    integer, intent(in) :: nume_inst
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hat_valinc(*)
    character(len=19), intent(in) :: sddyna
    integer, intent(in) :: list_func_acti(*)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! XFEM method - Initializations
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! IO  ds_contact       : datastructure for contact management
! In  nume_inst        : index of current step time
! In  hat_valinc       : hat variable for algorithm fields
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddyna           : datastructure for dynamic
! In  list_func_acti   : list of active functionnalities
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_dyna, l_step_first
    character(len=19) :: sdcont_depgeo, sdcont_deplam
    character(len=19) :: disp_prev
    character(len=19) :: xseuco, xseucp
    character(len=19) :: xindco, xmemco, xindcp, xmemcp, xcohes, xcohep
!
! --------------------------------------------------------------------------------------------------
!
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')
    ASSERT(.not. l_dyna)
!
! - Using *_INIT options (like SEUIL_INIT)
!
    l_step_first = nume_inst .eq. 1
!
! - Get field names in hat-variables
!
    call nmchex(hat_valinc, 'VALINC', 'DEPMOI', disp_prev)
!
! - Management of status for time cut
!
    xindco = ds_contact%sdcont_solv(1:14)//'.XFIN'
    xmemco = ds_contact%sdcont_solv(1:14)//'.XMEM'
    xindcp = ds_contact%sdcont_solv(1:14)//'.XFIP'
    xmemcp = ds_contact%sdcont_solv(1:14)//'.XMEP'
    xseuco = ds_contact%sdcont_solv(1:14)//'.XFSE'
    xseucp = ds_contact%sdcont_solv(1:14)//'.XFSP'
    xcohes = ds_contact%sdcont_solv(1:14)//'.XCOH'
    xcohep = ds_contact%sdcont_solv(1:14)//'.XCOP'
    call copisd('CHAMP_GD', 'V', xindcp, xindco)
    call copisd('CHAMP_GD', 'V', xmemcp, xmemco)
    call copisd('CHAMP_GD', 'V', xseucp, xseuco)
    call copisd('CHAMP_GD', 'V', xcohep, xcohes)
!
! - Save displacements for geometric loop
!
    sdcont_depgeo = ds_contact%sdcont_solv(1:14)//'.DEPG'
    call copisd('CHAMP_GD', 'V', disp_prev, sdcont_depgeo)
!
! - Save displacements for friction loop
!
    sdcont_deplam = ds_contact%sdcont_solv(1:14)//'.DEPF'
    call copisd('CHAMP_GD', 'V', disp_prev, sdcont_deplam)
!
! - Geometric loop counter initialization
!
    call mmbouc(ds_contact, 'Geom', 'Init_Counter')
!
! - First geometric loop counter
!
    call mmbouc(ds_contact, 'Geom', 'Incr_Counter')
!
! - Create fields
!
    call xmelem(mesh, model, ds_contact, list_func_acti)
!
end subroutine
