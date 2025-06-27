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
subroutine te0167(option, nomte)
!
    use Behaviour_module, only: behaviourOption
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
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUI_ABSO, AXIS_FLUI_ABSO
!
! Options: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: nx, ny, d(6), us(3), e(3, 3)
    real(kind=8) :: poids
    real(kind=8) :: rho, alpha, r_impe, rhon, q_alpha, q_c
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jv_geom, jv_mate, jv_matr
    integer(kind=8) :: jvVect, jvDisp, jvDispm, jvDispp
    integer(kind=8) :: ndim, nno, ndi, ipg, npg, n1, n2, nn
    integer(kind=8) :: ldec
    integer(kind=8) :: i, ij, j
    integer(kind=8) :: j_mater, iret, codret
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: FEForm, rela_comp
    aster_logical :: l_axis
    real(kind=8) :: r
    aster_logical :: lVect, lMatr, lVari, lSigm
!
! --------------------------------------------------------------------------------------------------
!
    lVect = ASTER_FALSE
    lMatr = ASTER_FALSE
    lVari = ASTER_FALSE
    lSigm = ASTER_FALSE
!
! - Check behaviour
!
    if (option(1:9) .eq. 'FULL_MECA' .or. &
        option .eq. 'RAPH_MECA' .or. &
        option .eq. 'RIGI_MECA_TANG') then
        call jevech('PCOMPOR', 'L', vk16=compor)
! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)
        rela_comp = compor(RELA_NAME)
        if (rela_comp .ne. 'ELAS') then
            call utmess('F', 'FLUID1_1')
        end if
    end if

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
!   for u-p-phi : K = 0
!   for u-p     : K =  rho  /Zr =  1.d0/(r_impe*q_alpha)
!   for u-psi   : K = -rho^2/Zr = -rho/(r_impe*q_alpha)
    rhon = -rho
    q_alpha = (1.d0+alpha)/(1.d0-alpha)
    q_c = r_impe*q_alpha
!
! - Compute
!
    if (FEForm .eq. 'U_P_PHI') then
        goto 110
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then

        do i = 1, ndi
            d(i) = 0.d0
        end do

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

            if (FEForm .eq. 'U_P') then
                if (r_impe .le. r8prem()) then
                    goto 110
                else
                    do i = 1, nno
                        do j = 1, i
                            ij = (i-1)*i/2+j
                            d(ij) = d(ij)+poids* &
                                    1.d0/q_c* &
                                    zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                        end do
                    end do
                end if
            elseif (FEForm .eq. 'U_PSI') then
                if (r_impe .le. r8prem()) then
                    goto 110
                else
                    do i = 1, nno
                        do j = 1, i
                            ij = (i-1)*i/2+j
                            d(ij) = d(ij)+poids* &
                                    rhon/q_c* &
                                    zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)
                        end do
                    end do
                end if
            end if
        end do
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if

! --------- Compute matrix
110 continue

    if (option(1:9) .eq. 'RIGI_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        call jevech('PMATUUR', 'E', jv_matr)

        do i = 1, ndi
            zr(jv_matr+i-1) = 0.d0
        end do

        if (FEForm .eq. 'U_P_PHI') then
            goto 120
        elseif (FEForm .eq. 'U_P') then
            if (r_impe .le. r8prem()) then
                goto 120
            else
                do i = 1, ndi
                    zr(jv_matr+i-1) = d(i)
                end do
            end if
        elseif (FEForm .eq. 'U_PSI') then
            if (r_impe .le. r8prem()) then
                goto 120
            else
                do i = 1, ndi
                    zr(jv_matr+i-1) = d(i)
                end do
            end if
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if

    if (option .eq. 'RIGI_MECA_HYST') then
        call jevech('PMATUUC', 'E', jv_matr)

        do i = 1, ndi
            zr(jv_matr+i-1) = dcmplx(0.d0, 0.d0)
        end do

        if (FEForm .eq. 'U_P_PHI') then
            goto 120
        elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            if (r_impe .le. r8prem()) then
                goto 120
            else
                do i = 1, ndi
                    zc(jv_matr+i-1) = dcmplx(d(i), 0.d0)
                end do
            end if
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if

!
! - Save vector
!
    if (lVect .or. option .eq. 'FORC_NODA') then
        call jevech('PVECTUR', 'E', jvVect)
        if (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            if (lVect) then
                call jevech('PDEPLMR', 'L', jvDispm)
                call jevech('PDEPLPR', 'L', jvDispp)
                do i = 1, nno
                    zr(jvVect+i-1) = 0.d0
                    us(i) = zr(jvDispm+i-1)+zr(jvDispp+i-1)
                end do
            elseif (option .eq. "FORC_NODA") then
                call jevech('PDEPLAR', 'L', jvDisp)
                do i = 1, nno
                    zr(jvVect+i-1) = 0.d0
                    us(i) = zr(jvDisp+i-1)
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
            nn = 0
            do n1 = 1, nno
                do n2 = 1, n1
                    nn = nn+1
                    e(n1, n2) = d(nn)
                    e(n2, n1) = d(nn)
                end do
            end do
            do n1 = 1, nno
                do n2 = 1, nno
                    zr(jvVect+n1-1) = zr(jvVect+n1-1)+e(n1, n2)*us(n2)
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if
120 continue
end subroutine
