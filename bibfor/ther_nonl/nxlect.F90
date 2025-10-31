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
subroutine nxlect(result, model, &
                  ther_crit_i, ther_crit_r, &
                  ds_inout, ds_algopara, &
                  ds_algorom, ds_print, &
                  compor, &
                  mesh, l_dry)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nxdocc.h"
#include "asterfort/nxdocn.h"
#include "asterfort/nxdomt.h"
#include "asterfort/dismoi.h"
#include "asterfort/nonlinDSInOutRead.h"
#include "asterfort/nonlinDSPrintRead.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: result
    character(len=8), intent(in) :: model
    integer(kind=8), intent(inout) :: ther_crit_i(*)
    real(kind=8), intent(inout) :: ther_crit_r(*)
    type(NL_DS_InOut), intent(inout) :: ds_inout
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=24), intent(out) :: compor
    character(len=8), intent(out) :: mesh
    aster_logical, intent(in) :: l_dry
!
!
    character(len=3) :: answer
!
! --------------------------------------------------------------------------------------------------
!
! Thermics - Init
!
! Read parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of datastructure for results
! In  model            : name of model
! IO  ther_crit_i      : criteria for algorithm (integer)
! IO  ther_crit_r      : criteria for algorithm (real)
! IO  ds_inout         : datastructure for input/output management
! IO  ds_algopara      : datastructure for algorithm parameters
! IO  ds_algorom       : datastructure for ROM parameters
! IO  ds_print         : datastructure for printing parameters
! Out compor           : name of <CARTE> COMPOR
! Out mesh             : name of mesh
! In  l_dry            : .true. if drying
!
! --------------------------------------------------------------------------------------------------
!
    compor = ' '
    mesh = ' '
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
! - Check model vs command
!
    if (l_dry) then
        call dismoi('EXI_NON_SECH', model, 'MODELE', repk=answer)
        if (answer .eq. 'OUI') call utmess('F', 'THERNONLINE4_4')
    else
        call dismoi('EXI_SECH', model, 'MODELE', repk=answer)
        if (answer .eq. 'OUI') call utmess('F', 'THERNONLINE4_3')
    end if
!
! - Create comportment <CARTE>
!
    compor = '&&NXDOCC.COMPOR'
    call nxdocc(model, compor)
!
! - Read parameters for algorithm management
!
    call nxdomt(ds_algopara, ds_algorom)
!
! - Read convergence criteria
!
    call nxdocn(ther_crit_i, ther_crit_r)
!
! - Read parameters for input/output management
!
    if (l_dry) then
        call nonlinDSInOutRead('SECH', result, ds_inout)
    else
        call nonlinDSInOutRead('THER', result, ds_inout)
    end if
!
! - Read parameters for printing
!
    call nonlinDSPrintRead(ds_print)
!
end subroutine
