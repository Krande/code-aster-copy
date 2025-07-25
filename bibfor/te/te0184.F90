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
subroutine te0184(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/teattr.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
#include "asterfort/utmess.h"
#include "asterfort/getFluidPara.h"
#include "asterc/r8prem.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUI_ABSO, AXIS_FLUI_ABSO
!
! Options: MASS_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nx, ny
    real(kind=8) :: poids
    real(kind=8) :: rho, alpha, r_impe, rhon, q_alpha, q_c
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jv_geom, jv_mate, jv_matr
    integer(kind=8) :: ndim, nno, ndi, ipg, npg
    integer(kind=8) :: ldec
    integer(kind=8) :: i, ii, ij, j, jj
    integer(kind=8) :: j_mater, iret
    character(len=16) :: FEForm
    aster_logical :: l_axis
    real(kind=8) :: r
!
! --------------------------------------------------------------------------------------------------
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    l_axis = (lteatt('AXIS', 'OUI'))
!
! - Get parameters of element
!
    call elrefe_info(fami='RIGI', &
                     ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. 3)
    if (FEForm .eq. 'U_P_PHI') then
        ndi = nno*(2*nno+1)
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndi = nno*(nno+1)/2
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho, alpha_=alpha, r_=r_impe)
!
! - Conditions on fluid parameters
!
    if ((1.d0-alpha) .ge. r8prem()) then
        q_alpha = (1.d0+alpha)/(1.d0-alpha)
    else
        call utmess('F', 'FLUID1_5', sk='ALPHA')
    end if

    if ((abs(rho)-rho) .gt. r8prem()) then
        call utmess('F', 'FLUID1_7', sk='RHO')
    end if

    if ((abs(r_impe)-r_impe) .gt. r8prem()) then
        call utmess('F', 'FLUID1_7', sk='R')
    end if
!
! - Output field
!
!   for u-p-phi : M = -rho^2/Zr = -rho/(r_impe*q_alpha)
!   for u-p     : M =  0
!   for u-psi   : M =  0
    rhon = -rho
    q_alpha = (1.d0+alpha)/(1.d0-alpha)
    q_c = r_impe*q_alpha

    call jevech('PMATUUR', 'E', jv_matr)

    do i = 1, ndi
        zr(jv_matr+i-1) = 0.d0
    end do
!
! - Compute
!
    do ipg = 1, npg
        ldec = (ipg-1)*nno
! --------- Compute normal
        nx = 0.0d0
        ny = 0.0d0
        call vff2dn(ndim, nno, ipg, ipoids, idfde, &
                    zr(jv_geom), nx, ny, poids)
! --------- Radius for axisymmetric model
        if (l_axis) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(jv_geom+2*(i-1))*zr(ivf+ldec+i-1)
            end do
            poids = poids*r
        end if
! --------- Compute matrix
        if (FEForm .eq. 'U_P_PHI') then
            if (r_impe .le. r8prem()) then
                goto 120
            else
                do i = 1, nno
                    do j = 1, i
                        ii = 2*i
                        jj = 2*j
                        ij = (ii-1)*ii/2+jj
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                           poids* &
                                           rhon/q_c* &
                                           zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                    end do
                end do
            end if
        elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            goto 120
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if

120     continue

    end do
end subroutine
