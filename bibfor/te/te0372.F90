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
subroutine te0372(option, nomte)
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
#include "asterfort/utmess.h"
#include "asterc/r8prem.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUIDE, AXIS_FLUIDE (boundary)
!
! Options: ONDE_FLUI
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_geom, jv_mate, jv_onde, jv_matr
    real(kind=8) :: nx, ny
    real(kind=8) :: poids, rho, celer
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: ndi, nno, npg, ndim
    aster_logical :: l_axis
    integer(kind=8) :: i, ii, ij, j, jj, ipg, ldec
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
    if (FEForm .eq. 'U_P_PHI') then
        ndi = nno*(2*nno+1)
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndi = nno*(nno+1)/2
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho, cele_r_=celer)
!
! - Output field
!
    call jevech('PMATUUR', 'E', jv_matr)
    do i = 1, ndi
        zr(jv_matr+i-1) = 0.d0
    end do
!
! - Compute
!
    if (abs(zr(jv_onde)) .gt. r8prem()) then
        do ipg = 1, npg
            ldec = (ipg-1)*nno
! --------- Compute normal
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
            if (FEForm .eq. 'U_P_PHI') then
                do i = 1, nno
                    do j = 1, i
                        ii = 2*i
                        jj = 2*j
                        ij = (ii-1)*ii/2+jj
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)- &
                                           poids*zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)*rho/celer
                    end do
                end do
            elseif (FEForm .eq. 'U_P') then
                do i = 1, nno
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                           poids*zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)/celer
                    end do
                end do
            elseif (FEForm .eq. 'U_PSI') then
                do i = 1, nno
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)- &
                                           poids*zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)*rho/celer
                    end do
                end do
            else
                call utmess('F', 'FLUID1_2', sk=FEForm)
            end if
        end do
    end if
!
end subroutine
