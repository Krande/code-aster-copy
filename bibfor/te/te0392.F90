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

subroutine te0392(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/caatdb.h"
#include "asterfort/cast3d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elraga.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/invjac.h"
#include "asterfort/jevech.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/nbsigm.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
!
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_SI
! Option: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: idfde2, igau, imate, imatuu, ipoid2
    integer(kind=8) :: nbsig, nno, npg1
    real(kind=8) :: jacgau
    real(kind=8) :: angl_naut(3), instan
    integer(kind=8) :: igeom, ipoids, ivf, idfde
!
    aster_logical :: calbn
    integer(kind=8) :: i, ino, j, k, proj, nbpg2, ipg, ispg
    integer(kind=8) :: ndim, nnos, kp
    real(kind=8) :: d(6, 6), s
    real(kind=8) :: poipg2(8), b(6, 81), b0(6, 3, 8)
    real(kind=8) :: jac, invja(3, 3), bi(3, 8), hx(3, 4)
    real(kind=8) :: gam(4, 8), coopg2(24), h(8, 4), dh(4, 24)
    real(kind=8) :: bn(6, 3, 8)
    real(kind=8) :: dfdx(8), dfdy(8), dfdz(8)
    real(kind=8) :: nu, nub, nu12
    character(len=16) :: elas_keyword
    integer(kind=8) :: elas_id
    data h/1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0,&
     &        1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0,&
     &        1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0,&
     &       -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0, 1.d0, -1.d0/
!
! --------------------------------------------------------------------------------------------------
!
!
! - Finite element informations
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Initializations
!
    instan = r8vide()
    b(:, :) = 0.d0
!
! - Number of stress components
!
    nbsig = nbsigm()
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)
!
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
    call get_elas_id(zi(imate), elas_id, elas_keyword)
!
! - Orthotropic parameters
!
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
    call jevech('PMATUUR', 'E', imatuu)
    do i = 1, 300
        zr(imatuu-1+i) = 0.0d0
    end do
!
! - Compute [Bi] (mean value for derivate of shape functions)
!
    do ipg = 1, npg1
        call dfdm3d(nno, ipg, ipoids, idfde, zr(igeom), &
                    jac, dfdx, dfdy, dfdz)
        do ino = 1, nno
            bi(1, ino) = dfdx(ino)
            bi(2, ino) = dfdy(ino)
            bi(3, ino) = dfdz(ino)
        end do
    end do
!
    do igau = 1, npg1
!
! ----- Compute matrix [B]: displacement -> strain (first order)
!
        call dfdm3d(nno, igau, ipoids, idfde, zr(igeom), &
                    jacgau, dfdx, dfdy, dfdz)
!
! ----- Modify matrix [B] for underintegrated elements
!
        do i = 1, 8
            j = 3*(i-1)+1
            b(1, j) = bi(1, i)
            b(2, j+1) = bi(2, i)
            b(3, j+2) = bi(3, i)
            b(4, j) = bi(2, i)
            b(4, j+1) = bi(1, i)
            b(5, j) = bi(3, i)
            b(5, j+2) = bi(1, i)
            b(6, j+1) = bi(3, i)
            b(6, j+2) = bi(2, i)
        end do
        do i = 1, nno
            do j = 1, 3
                do k = 1, 6
                    b0(k, j, i) = b(k, (i-1)*3+j)
                end do
            end do
        end do
!
! ----- Compute Hooke matrix [D]
!
        call dmatmc('RIGI', zi(imate), instan, '+', igau, &
                    1, angl_naut, nbsig, d)
!
! ----- Compute "center" rigidity matrix [KC]
!
        call caatdb(nno, b0, d, b0, jacgau, &
                    zr(imatuu))
!
    end do
!
! - Gamma ratio
!
    do i = 1, 4
        do k = 1, 3
            hx(k, i) = 0.d0
            do j = 1, nno
                hx(k, i) = hx(k, i)+h(j, i)*zr(igeom-1+3*(j-1)+k)
            end do
        end do
    end do
!
    do i = 1, 4
        do j = 1, nno
            s = 0.d0
            do k = 1, 3
                s = s+hx(k, i)*bi(k, j)
            end do
            gam(i, j) = 0.125d0*(h(j, i)-s)
        end do
    end do
!
! - Poisson ration for ASQBI
!
    ipg = 1
    ispg = 1
    call get_elas_para('RIGI', zi(imate), '+', ipg, ispg, &
                       elas_id, elas_keyword, &
                       nu_=nu, nu12_=nu12)
    if (elas_id .eq. 1) then
        nub = nu/(1.d0-nu)
    else
        nub = nu12/(1.d0-nu12)
    end if
!
! - Projection type
!           0 AUCUNE
!           1 ADS
!           2 ASBQI
!
    proj = 2
    calbn = .false.
!
! - Finite element informations for underintegrated element
!
    call elraga('HE8', 'FPG8    ', ndim, nbpg2, coopg2, &
                poipg2)
    call elrefe_info(elrefe='HE8', fami='MASS', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=nbpg2, jpoids=ipoid2, jdfde=idfde2)
!
! - Compute corrected stabilization matrix [K_STAB]
!
    do ipg = 1, nbpg2
        kp = 3*(ipg-1)
        call invjac(nno, ipg, ipoid2, idfde2, zr(igeom), &
                    invja, jac)
        do i = 1, 3
            dh(1, kp+i) = coopg2(3*ipg-1)*invja(i, 3)+coopg2(3*ipg)*invja(i, 2)
        end do
        do i = 1, 3
            dh(2, kp+i) = coopg2(3*ipg-2)*invja(i, 3)+coopg2(3*ipg)*invja(i, 1)
        end do
        do i = 1, 3
            dh(3, kp+i) = coopg2(3*ipg-2)*invja(i, 2)+coopg2(3*ipg-1)*invja(i, 1)
        end do
        do i = 1, 3
            dh(4, kp+i) = coopg2(3*ipg-2)*coopg2(3*ipg-1)*invja(i, 3)+coopg2(3*ipg-1)*coop&
                         &g2(3*ipg)*invja(i, 1)+coopg2(3*ipg-2)*coopg2(3*ipg)*invja(i, 2)
        end do
        call cast3d(proj, gam, dh, b0, nno, &
                    ipg, nub, nu, d, calbn, &
                    bn, jac, zr(imatuu))
    end do
!
end subroutine
