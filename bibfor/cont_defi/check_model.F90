! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine check_model(mesh, model, cont_form)
!
implicit none
!
#include "asterfort/exixfe.h"
#include "asterfort/exipat.h"
#include "asterfort/dismoi.h"
#include "asterfort/utmess.h"
!
#ifdef _USE_MPI
#include "mpif.h"
#include "asterf_mpi.h"
#endif
!
!
character(len=8), intent(in) :: mesh, model    
integer, intent(in) :: cont_form
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Check if exist XFEM in model or PATCH in mesh
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  cont_form        : formulation of contact
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iret
!
! --------------------------------------------------------------------------------------------------

! - Check if exist XFEM in model
!
    if (cont_form .eq. 3) then
        call exixfe(model, iret)
        if (iret .eq. 0) then
            call utmess('F', 'XFEM2_8', sk=model)
        endif
    endif
!
! - Check if exist PATCH in mesh (LAC method)
    if (cont_form .eq. 5) then
        call exipat(mesh, iret)
        if (iret .eq. 0) then
            call utmess('F', 'CONTACT4_2', sk=mesh)
        endif
    endif  
end subroutine
 
