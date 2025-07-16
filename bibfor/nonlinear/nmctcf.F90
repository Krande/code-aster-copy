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
subroutine nmctcf(mesh, sderro, hval_incr, ds_print, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/cfdisi.h"
#include "asterfort/copisd.h"
#include "asterfort/infdbg.h"
#include "asterfort/mmbouc.h"
#include "asterfort/mmmcri_frot.h"
#include "asterfort/mmreas.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmimck.h"
#include "asterfort/nmimcr.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: hval_incr(*)
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algo
!
! Friction loop management - Management
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  sderro           : datastructure for errors during algorithm
! In  hval_incr        : hat-variable for incremental values fields
! IO  ds_print         : datastructure for printing parameters
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    aster_logical :: loop_fric_error
    integer(kind=8) :: iter_fric_maxi
    integer(kind=8) :: loop_fric_count
    character(len=19) :: disp_curr, loop_fric_disp
    character(len=16) :: loop_fric_node
    real(kind=8) :: loop_fric_vale
    aster_logical :: loop_fric_conv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> MISE A JOUR DU SEUIL DE TRESCA'
    end if
!
! - Initializations
!
    loop_fric_conv = .false.
    loop_fric_node = ' '
    loop_fric_vale = r8vide()
    call mmbouc(ds_contact, 'Fric', 'Set_NoError')
!
! - Get fields
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)

! - Get friction loop parameters
    loop_fric_disp = ds_contact%sdcont_solv(1:14)//'.DEPF'
    iter_fric_maxi = cfdisi(ds_contact%sdcont_defi, 'ITER_FROT_MAXI')
!
! - Update triggers
!
    call mmreas(mesh, ds_contact, hval_incr)

!
! - Compute friction criterion
!
    call mmmcri_frot(mesh, loop_fric_disp, disp_curr, ds_contact)
!
! - Get final loop state
!
    call mmbouc(ds_contact, 'Fric', 'Is_Convergence', loop_state_=loop_fric_conv)
    call mmbouc(ds_contact, 'Fric', 'Get_Locus', loop_locus_=loop_fric_node)
    call mmbouc(ds_contact, 'Fric', 'Get_Vale', loop_vale_=loop_fric_vale)
    call mmbouc(ds_contact, 'Fric', 'Read_Counter', loop_fric_count)
!
! - Too many iterations ?
!
    if ((.not. loop_fric_conv) .and. (loop_fric_count .eq. iter_fric_maxi)) then
        call mmbouc(ds_contact, 'Fric', 'Set_Error')
    end if
    call mmbouc(ds_contact, 'Fric', 'Is_Error', loop_state_=loop_fric_error)
!
! - Save events
!
    call nmcrel(sderro, 'ERRE_CTCF', loop_fric_error)
    if (loop_fric_conv) then
        call nmcrel(sderro, 'DIVE_FIXF', .false._1)
    else
        call nmcrel(sderro, 'DIVE_FIXF', .true._1)
    end if
!
! - Set values in convergence table for contact geoemtry informations
!
    call nmimck(ds_print, 'BOUC_NOEU', loop_fric_node, .true._1)
    call nmimcr(ds_print, 'BOUC_VALE', loop_fric_vale, .true._1)
!
! - Update reference displacement for friction loop
!
    if (.not. loop_fric_conv) then
        call copisd('CHAMP_GD', 'V', disp_curr, loop_fric_disp)
    end if
!
end subroutine
