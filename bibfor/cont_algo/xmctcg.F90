! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine xmctcg(model, mesh, ds_contact, ds_measure)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/xappar.h"
#include "asterfort/xreacg.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! XFEM method - Geometric loop: geometric actualisation and pairing
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm,*) '<CONTACT> . Geometric actualisation and pairing'
    endif
!
! - Geometric loop: new geometric iteration for statistics
!
    call nmrinc(ds_measure, 'Cont_Geom')
!
! - Geometric loop: begin timer
!
    call nmtime(ds_measure, 'Init'  , 'Cont_Geom')
    call nmtime(ds_measure, 'Launch', 'Cont_Geom')
!
! - Geometric actualisation
!
    call xreacg(model, ds_contact)
!
! - Pairing
!
    call xappar(mesh, model, ds_contact)
!
! - Geometric loop: end timer
!
    call nmtime(ds_measure, 'Stop', 'Cont_Geom')
!
end subroutine
