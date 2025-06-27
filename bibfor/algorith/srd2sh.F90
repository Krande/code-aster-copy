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

subroutine srd2sh(nmat, materf, varh, dhds, devsig, rcos3t, d2shds)

!

!!!
!!! MODELE LKR : CALCUL DE DERIVEE 2NDE DE SII*H PAR RAPPORT A SIGMA
!!!

! ===================================================================================
! IN  : NMAT           : DIMENSION TABLE DES PARAMETRES MATERIAU
!     : MATERF(NMAT,2) : PARAMETRES MATERIAU A T+DT
!     : VARH           : VECTEUR CONTENANT H0E,H0C ET HTHETA
!     : DHDS(6)        : DERIVEE DE HTHETA PAR RAPPORT A SIGMA
!     : DEVSIG(6)      : DEIATEUR DES CONTRAINTES
!     : RCOS3T         : COS(3THETA) = SQRT(54)*DET(DEVISG)/SII**3
! OUT : D2SHDS(6,6)    :  DERIVEE 2NDE SII*H PAR RAPPORT A SIGMA (NDT X NDT)
!     : IRET           :  CODE RETOUR
! ===================================================================================

    implicit none

#include "asterfort/lcprte.h"
#include "asterfort/srd2hs.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nmat
    real(kind=8) :: materf(nmat, 2), varh(2), d2shds(6, 6), dhds(6), devsig(6), rcos3t

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: ndi, ndt, i, j
    real(kind=8) :: sii, dikdjl(6, 6), dijdkl(6, 6)
    real(kind=8) :: dsiids(6), dsdsig(6, 6), mat1(6, 6), d2hds2(6, 6)
    real(kind=8) :: mat2(6, 6), mat3(6, 6), dhtds(6), mat4(6, 6), mat5(6, 6)
    real(kind=8) :: d2hdsi(6, 6)
    common/tdim/ndt, ndi

    !!!
    !!! Construction de sii
    !!!

    sii = sqrt(dot_product(devsig(1:ndt), devsig(1:ndt)))

    !!!
    !!! initialisation matrice d_ik x d_jl
    !!!

    dikdjl(:, :) = 0.d0

    do i = 1, ndt
        dikdjl(i, i) = 1.d0
    end do

    !!!
    !!! initialisation matrice d_ij x d_kl
    !!!

    dijdkl(:, :) = 0.d0

    do i = 1, ndi
        do j = 1, ndi
            dijdkl(i, j) = 1.d0/3.d0
        end do
    end do

    !!!
    !!! Calcul de d(sii)/d(sigma)
    !!!

    do i = 1, ndt
        dsiids(i) = 0.d0
        do j = 1, ndt
            dsdsig(j, i) = dikdjl(j, i)-dijdkl(j, i)
            dsiids(i) = dsiids(i)+devsig(j)*dsdsig(j, i)/sii
        end do
    end do

    !!! Calcul de dh/ds*dsii/dsig = mat1
    call lcprte(dhds, dsiids, mat1)

    !!! Calcul de d2h/ds2
    call srd2hs(nmat, materf, devsig, sii, rcos3t, d2hds2)

    !!! Calcul de d2h/ds2
    d2hdsi(1:ndt, 1:ndt) = matmul(d2hds2(1:ndt, 1:ndt), dsdsig(1:ndt, 1:ndt))

    !!! Construction de sii*d2h/dsigma = mat2
    mat2(1:ndt, 1:ndt) = sii*d2hdsi(1:ndt, 1:ndt)

    !!! mat2 + mat1 = mat3
    mat3(1:ndt, 1:ndt) = mat1(1:ndt, 1:ndt)+mat2(1:ndt, 1:ndt)

    !!! mat2 = coefh*mat3
    mat2(1:ndt, 1:ndt) = mat3(1:ndt, 1:ndt)

    !!! Construction de dh/dsigma = (dh/ds)*(ds/dsigma)
    mat1(1:ndt, 1:ndt) = dsdsig(1:ndt, 1:ndt)

    do i = 1, ndt
        dhtds(i) = 0.d0
        do j = 1, ndt
            dhtds(i) = dhtds(i)+dhds(j)*mat1(j, i)/sii
        end do
    end do

    !!! Construction de mat1 = s x dh/dsigma
    call lcprte(devsig, dhtds, mat1)

    !!! Calcul de mat3 ) h/sii*(ds/dsig)
    mat3(1:ndt, 1:ndt) = (varh(2)/sii)*dsdsig(1:ndt, 1:ndt)

    !!! Calcul de mat5 = h*s*(dsii/ds)/sii**2
    call lcprte(devsig, dsiids, mat4)
    mat5(1:ndt, 1:ndt) = (varh(2)/sii**2.d0)*mat4(1:ndt, 1:ndt)

    !!! mat4 = mat2+mat1+mat3-mat5
    do i = 1, ndt
        do j = 1, ndt
            mat4(i, j) = mat2(i, j)+mat1(i, j)+mat3(i, j)-mat5(i, j)
        end do
    end do

    !!! d2sh/ds = ds/dsig * mat4
    d2shds(1:ndt, 1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), mat4(1:ndt, 1:ndt))

end subroutine
