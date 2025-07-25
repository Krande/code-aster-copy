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
subroutine thmCompVariElno(ds_thm)
!
    use Behaviour_type
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/thmGetElemRefe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/posthm.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/thmGetElemIntegration.h"
#include "asterfort/Behaviour_type.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! VARI_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'VARI_ELNO'
    character(len=8) :: elrefe, elref2
    integer(kind=8) :: jv_gano, jv_compo, jv_varielga, jv_varielno
    integer(kind=8) :: ncmp, nvim
    aster_logical :: l_axi, l_vf
    character(len=3) :: inte_type
    integer(kind=8) :: ndim
!
! --------------------------------------------------------------------------------------------------
!

! - Get model of finite element
    call thmGetElemModel(ds_thm, l_axi, l_vf, ndim)
!
! - Get type of integration
!
    call thmGetElemIntegration(l_vf, inte_type)
!
! - Get reference elements
!
    call thmGetElemRefe(l_vf, elrefe, elref2)
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, jgano=jv_gano)
!
! - Parameters from behaviour
!
    call jevech('PCOMPOR', 'L', jv_compo)
    read (zk16(jv_compo-1+NVAR), '(I16)') ncmp
    read (zk16(jv_compo-1+MECA_NVAR), '(I16)') nvim
!
! - Compute
!
    call jevech('PVARIGR', 'L', jv_varielga)
    call jevech('PVARINR', 'E', jv_varielno)
    call posthm(option, inte_type, jv_gano, ncmp, nvim, &
                zr(jv_varielga), zr(jv_varielno))
!
end subroutine
