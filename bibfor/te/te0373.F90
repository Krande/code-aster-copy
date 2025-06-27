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
subroutine te0373(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/vff2dn.h"
#include "asterfort/teattr.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUIDE, AXIS_FLUIDE (boundary)
!
! Options: CHAR_MECA_ONDE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_geom, jv_mate, jv_onde, jv_vect
    real(kind=8) :: nx, ny
    real(kind=8) :: poids, celer
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: nno, npg, ndim, ndofbynode
    integer(kind=8) :: i, ii, ipg, ldec
    aster_logical :: l_axis
    real(kind=8) :: r
    integer(kind=8) :: j_mater, iret
    character(len=16) :: FEForm
!
! --------------------------------------------------------------------------------------------------
!

!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
    call jevech('PONDECR', 'L', jv_onde)
!
! - Get element parameters
!
    l_axis = (lteatt('AXIS', 'OUI'))
    call teattr('S', 'FORMULATION', FEForm, iret)
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. 3)
    if (FEForm .eq. 'U_P_PHI') then
        ndofbynode = 2
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndofbynode = 1
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, cele_r_=celer)
!
! - Output field
!
    call jevech('PVECTUR', 'E', jv_vect)
    do i = 1, ndofbynode*nno
        zr(jv_vect+i-1) = 0.d0
    end do
!
! - Loop on Gauss points
!
    do ipg = 1, npg
        ldec = (ipg-1)*nno
! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        call vff2dn(ndim, nno, ipg, ipoids, idfde, &
                    zr(jv_geom), nx, ny, poids)
        if (l_axis) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(jv_geom+2*(i-1))*zr(ivf+ldec+i-1)
            end do
            poids = poids*r
        end if
! ----- Compute vector
        do i = 1, nno
            ii = ndofbynode*i
            zr(jv_vect+ii-1) = zr(jv_vect+ii-1)+ &
                               poids*zr(jv_onde+ipg-1)*zr(ivf+ldec+i-1)/celer
        end do
    end do
!
end subroutine
