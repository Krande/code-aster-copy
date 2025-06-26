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
subroutine te0253(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: AXIS_FLUIDE, 2D_FLUIDE
!
! Options: RIGI_MECA/FORC_NODA/FULL_MECA/RAPH_MECA/RIGI_MECA_HYST/RIGI_MECA_TANG
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l
    integer(kind=8) :: n1, n2
    integer(kind=8) :: nn, nno2, nt1, nt2
    integer(kind=8) :: ipg, ij, ik, ijkl, codret
    integer(kind=8) :: jvDispm, jvDispp, jvDisp
    integer(kind=8) :: jv_geom, jv_mate
    integer(kind=8) :: jvVect, jv_codret, jv_matr
    character(len=16) :: rela_comp
    real(kind=8) :: a(2, 2, 9, 9), mmat(9, 9)
    real(kind=8) :: b(18, 18), e(9, 9), ul(18), us(9), c(171), d(45)
    real(kind=8) :: dfdx(9), dfdy(9)
    real(kind=8) :: poids, rho, celer
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: nno, npg
    integer(kind=8) :: j_mater, iret
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: FEForm
    aster_logical :: l_axis
    real(kind=8) :: r
    aster_logical :: lVect, lMatr, lVari, lSigm
!
! --------------------------------------------------------------------------------------------------
!
    a = 0.d0
    mmat = 0.d0
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
!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
!
! - Get element parameters
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    l_axis = (lteatt('AXIS', 'OUI'))
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. 9)
    nno2 = nno*2
    nt2 = nno*(nno2+1)
    nt1 = nno*(nno+1)/2
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho, cele_r_=celer)
!
! - Loop on Gauss points
!
    do ipg = 1, npg
        k = (ipg-1)*nno
        call dfdm2d(nno, ipg, ipoids, idfde, zr(jv_geom), &
                    poids, dfdx, dfdy)
        if (l_axis) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(jv_geom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
        if (FEForm .eq. 'U_P_PHI') then
            do i = 1, nno
                do j = 1, i
                    if (abs(celer) .le. r8prem() .or. abs(rho) .le. r8prem()) then
                        a(1, 1, i, j) = 0.d0
                    else
                        a(1, 1, i, j) = a(1, 1, i, j)+ &
                                        poids*zr(ivf+k+i-1)*zr(ivf+k+j-1)/rho/celer**2.d0
                    end if
                end do
            end do
        elseif (FEForm .eq. 'U_P') then
            do j = 1, nno
                do i = 1, j
                    if (abs(celer) .le. r8prem() .or. abs(rho) .le. r8prem()) then
                        mmat(i, j) = 0.d0
                    else
                        mmat(i, j) = mmat(i, j)+ &
                                     poids*(dfdx(i)*dfdx(j)+dfdy(i)*dfdy(j))
                    end if
                end do
            end do
        elseif (FEForm .eq. 'U_PSI') then
            do j = 1, nno
                do i = 1, j
                    if (abs(celer) .le. r8prem() .or. abs(rho) .le. r8prem()) then
                        mmat(i, j) = 0.d0
                    else
                        mmat(i, j) = mmat(i, j)-rho* &
                                     poids*(dfdx(i)*dfdx(j)+dfdy(i)*dfdy(j))
                    end if
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end do
!
! - Compute result
!
    if (FEForm .eq. 'U_P_PHI') then
        do k = 1, 2
            do l = 1, 2
                do i = 1, nno
                    ik = ((2*i+k-3)*(2*i+k-2))/2
                    do j = 1, i
                        ijkl = ik+2*(j-1)+l
                        c(ijkl) = a(k, l, i, j)
                    end do
                end do
            end do
        end do
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        do j = 1, nno
            do i = 1, j
                ij = (j-1)*j/2+i
                d(ij) = mmat(i, j)
            end do
        end do
    end if
!
! - Save matrix
!
    if (option .eq. 'RIGI_MECA_HYST') then
        call jevech('PMATUUC', 'E', jv_matr)
        if (FEForm .eq. 'U_P_PHI') then
            do i = 1, nt2
                zc(jv_matr+i-1) = dcmplx(c(i), 0.d0)
            end do
        elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            do j = 1, nno
                do i = 1, j
                    ij = (j-1)*j/2+i
                    zc(jv_matr+ij-1) = dcmplx(mmat(i, j), 0.d0)
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    elseif (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA') then
        call jevech('PMATUUR', 'E', jv_matr)
        if (FEForm .eq. 'U_P_PHI') then
            do i = 1, nt2
                zr(jv_matr+i-1) = c(i)
            end do
        elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            do i = 1, nt1
                zr(jv_matr+i-1) = d(i)
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
        ! elseif (option(1:9) .eq. 'RIGI_MECA') then
        !     call jevech('PMATUUR', 'E', jv_matr)
        !     if (FEForm .eq. 'U_P_PHI') then
        !         do i = 1, nt2
        !             zr(jv_matr+i-1) = c(i)
        !         end do
        !     elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        !         do j = 1, nno
        !             do i = 1, j
        !                 ij = (j-1)*j/2+i
        !                 zr(jv_matr+ij-1) = mmat(i, j)
        !             end do
        !         end do
        !     else
        !         call utmess('F', 'FLUID1_2', sk=FEForm)
        !     end if
    end if
!
! - Save vector
!
    if (lVect .or. option .eq. 'FORC_NODA') then
        call jevech('PVECTUR', 'E', jvVect)
        if (FEForm .eq. 'U_P_PHI') then
            if (lVect) then
                call jevech('PDEPLMR', 'L', jvDispm)
                call jevech('PDEPLPR', 'L', jvDispp)
                do i = 1, nno2
                    zr(jvVect+i-1) = 0.d0
                    ul(i) = zr(jvDispm+i-1)+zr(jvDispp+i-1)
                end do
            elseif (option .eq. "FORC_NODA") then
                call jevech('PDEPLAR', 'L', jvDisp)
                do i = 1, nno2
                    zr(jvVect+i-1) = 0.d0
                    ul(i) = zr(jvDisp+i-1)
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
            nn = 0
            do n1 = 1, nno2
                do n2 = 1, n1
                    nn = nn+1
                    b(n1, n2) = c(nn)
                    b(n2, n1) = c(nn)
                end do
            end do
            do n1 = 1, nno2
                do n2 = 1, nno2
                    zr(jvVect+n1-1) = zr(jvVect+n1-1)+b(n1, n2)*ul(n2)
                end do
            end do
        elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
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
!
! - Save return code
!
    if (lSigm) then
        if (FEForm .eq. 'U_P_PHI' .or. FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
            call jevech('PCODRET', 'E', jv_codret)
            zi(jv_codret) = 0
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if
!
end subroutine
