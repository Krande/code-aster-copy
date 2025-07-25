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
subroutine diinit(mesh_, model_, ds_inout, mate, mateco, cara_elem, &
                  list_func_acti, sddyna, ds_conv, ds_algopara, solver, &
                  ds_contact, sddisc)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/getvid.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmcrar.h"
#include "asterfort/nmcrli.h"
#include "asterfort/nmcrsu.h"
#include "asterfort/ndxcfl.h"
!
    character(len=*), intent(in) :: mesh_
    character(len=*), intent(in) :: model_
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: cara_elem
    character(len=24), intent(in) :: mate, mateco
    type(NL_DS_Conv), intent(in) :: ds_conv
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: solver
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: list_func_acti(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Datastructures
!
! Time discretization and storing datastructures
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  mate             : name of material characteristics (field)
! In  cara_elem        : name of elementary characteristics (field)
! In  sddyna           : name of dynamic parameters
! In  ds_conv          : datastructure for convergence management
! In  ds_algo          : datastructure for algorithm management
! In  list_func_acti   : active functionnalities vector (see nmfonc)
! In  ds_contact       : datastructure for contact management
! In  solver           : name of solver parameters
! In  sddisc           : datastructure for time discretization
! In  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_expl, l_implex
    character(len=19) :: listInst
    character(len=8) :: model, mesh, result
!
! --------------------------------------------------------------------------------------------------
!
    call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=listInst)
    model = model_
    mesh = mesh_

! - Get parameters
    result = ds_inout%result

! - Active functionnalities
    l_expl = ndynlo(sddyna, 'EXPLICITE')
    l_implex = isfonc(list_func_acti, 'IMPLEX')

! - Create time discretization datastructure
    call nmcrli(listInst, sddisc)

! - Courant condition
    if (l_expl) then
        call ndxcfl(mate, mateco, cara_elem, sddyna, sddisc)
    end if

! - Create storing datastructure
    call nmcrar(result, sddisc, list_func_acti)

! - Automatic management of time stepping
    call nmcrsu(sddisc, listInst, ds_conv, ds_algopara, l_implex, &
                solver, ds_contact)
!
end subroutine
