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
subroutine te0258(option, nomte)
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
! Options: AMOR_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nx, ny, jac
    real(kind=8) :: poids
    real(kind=8) :: rho, rhon, cele_r, alpha, q_alpha, q_c, coef_ordre, onde_flui
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jv_geom, jv_mate, jv_matr, jv_amor
    integer(kind=8) :: ndim, nno, ndi, ipg, npg
    integer(kind=8) :: ldec
    integer(kind=8) :: i, ij, j
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
        ndi = 4*nno*nno
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
    call jevech('PAMORFL', 'L', jv_amor)
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho, cele_r_=cele_r, alpha_=alpha)
!
! - Conditions on fluid parameters
!
    if ((1.d0-alpha) .ge. r8prem()) then
        q_alpha = (1.d0+alpha)/(1.d0-alpha)
    else
        call utmess('F', 'FLUID1_5', sk='ALPHA')
    end if

    if (rho .le. r8prem()) then
        call utmess('F', 'FLUID1_6', sk='RHO')
    end if

    if (cele_r .le. r8prem()) then
        call utmess('F', 'FLUID1_6', sk='CELE_R')
    end if

    if (FEForm .eq. 'U_P_PHI') then
        if (zi(jv_amor-1+1) .eq. 1) then
            coef_ordre = 1.d0
        else
            coef_ordre = 0.d0
        end if
    else
        coef_ordre = 1.d0
    end if

    if (zi(jv_amor-1+2) .eq. 1) then
        onde_flui = +1.d0
    else
        onde_flui = -1.d0
    end if
!
! - Output field
!
!   for u-p-phi : Q =  rho  /Zc =  1.d0/(cele_r*q_alpha)
!   for u-p     : Q =  rho  /Zc =  1.d0/(cele_r*q_alpha)
!   for u-psi   : Q = -rho^2/Zc = -rho/(cele_r*q_alpha)
    rhon = -rho
    q_alpha = (1.d0+alpha)/(1.d0-alpha)
    q_c = cele_r*q_alpha

    if (FEForm .eq. 'U_P_PHI') then
        call jevech('PMATUNS', 'E', jv_matr)
    elseif (FEForm .eq. 'U_P') then
        call jevech('PMATUUR', 'E', jv_matr)
    elseif (FEForm .eq. 'U_PSI') then
        call jevech('PMATUUR', 'E', jv_matr)
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if

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
! --------- Compute jacobian
        jac = sqrt(nx*nx+ny*ny)
! --------- Compute matrix
        if (FEForm .eq. 'U_P_PHI') then
            do i = 1, nno
                do j = 1, nno
                    ij = 4*nno*(i-1)+2*j
                    zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                       jac*poids* &
                                       onde_flui*coef_ordre/q_c* &
                                       zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                end do
            end do
        elseif (FEForm .eq. 'U_P') then
            do i = 1, nno
                do j = 1, i
                    ij = (i-1)*i/2+j
                    zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                       jac*poids* &
                                       onde_flui*coef_ordre/q_c* &
                                       zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                end do
            end do
        elseif (FEForm .eq. 'U_PSI') then
            do i = 1, nno
                do j = 1, i
                    ij = (i-1)*i/2+j
                    zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                       jac*poids* &
                                       onde_flui*coef_ordre*rhon/q_c* &
                                       zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end do
end subroutine
