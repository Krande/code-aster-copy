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
subroutine te0420(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/teattr.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "asterfort/elrefe_info.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_FLUIDE, AXIS_FLUIDE, 2D_FLUIDE
!
! Options: PRME_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=16) :: FEForm
    integer(kind=8) :: nno, ino, ndofbynode
    integer(kind=8) :: jv_prme, jv_pres
!
! --------------------------------------------------------------------------------------------------
!
    nno = 0
!
! - Get element parameters
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    call elrefe_info(fami='RIGI', nno=nno)
    if (FEForm .eq. 'U_P_PHI') then
        ndofbynode = 2
    elseif (FEForm .eq. 'U_P') then
        ndofbynode = 1
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Input field
!
    call jevech('PDEPLAC', 'L', jv_pres)
!
! - Output field
!
    call jevech('PPRME_R', 'E', jv_prme)
    do ino = 1, nno
        zr(jv_prme+ino-1) = 20.d0*log10(abs(zc(jv_pres+ndofbynode*(ino-1)))/2.d-5)
    end do
!
end subroutine
