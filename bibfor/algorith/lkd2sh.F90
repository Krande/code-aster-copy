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
subroutine lkd2sh(nmat, materf, varh, dhds, devsig, &
                  rcos3t, d2shds, iret)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!     ------------------------------------------------------------------
!     CALCUL DE DERIVEE 2NDE DE SII*H PAR RAPPORT A SIGMA
!     IN  NMAT   : DIMENSION TABLE DES PARAMETRES MATERIAU
!         MATERF : PARAMETRES MATERIAU A T+DT
!         VARH   : VECTEUR CONTENANT H0E,H0C ET HTHETA
!         DHDS   : DERIVEE DE HTHETA PAR RAPPORT A SIGMA
!         DEVSIG : DEIATEUR DES CONTRAINTES
!         RCOS3T : COS(3THETA) = SQRT(54)*DET(DEVISG)/SII**3
!     OUT D2SHDS :  DERIVEE 2NDE SII*H PAR RAPPORT A SIGMA (NDT X NDT)
!         IRET   :  CODE RETOUR
!     ------------------------------------------------------------------
#include "asterfort/lcprte.h"
#include "asterfort/lkd2hs.h"
    integer(kind=8) :: iret, nmat
    real(kind=8) :: materf(nmat, 2), varh(3), d2shds(6, 6), dhds(6)
    real(kind=8) :: devsig(6), rcos3t
!
    integer(kind=8) :: ndi, ndt, i, j
    real(kind=8) :: h0ext, coefh, sii, un, zero, dikdjl(6, 6), dijdkl(6, 6)
    real(kind=8) :: trois, dsiids(6), dsdsig(6, 6), mat1(6, 6), d2hds2(6, 6)
    real(kind=8) :: mat2(6, 6), mat3(6, 6), dhtds(6), mat4(6, 6), mat5(6, 6)
    real(kind=8) :: d2hdsi(6, 6)
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(trois=3.0d0)
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!     ------------------------------------------------------------------
!
! --- RECUPERATION PROPRIETES MATERIAUX
    h0ext = materf(4, 2)
!
! --- COEFFICIENT (H0C-H0EXT)/(H0C-HOE)
    coefh = (varh(2)-h0ext)/(varh(2)-varh(1))
!
! --- CONSTRUCTION DE SII
    sii = norm2(devsig(1:ndt))
!
! --- INITIALISATION MATRICE D_IK X D_JL
    dikdjl(:, :) = zero
    do i = 1, ndt
        dikdjl(i, i) = un
    end do
!
! --- INITIALISATION MATRICE D_IJ X D_KL
    dijdkl(:, :) = zero
    do i = 1, ndi
        do j = 1, ndi
            dijdkl(i, j) = un/trois
        end do
    end do
!
! --- CALCUL DERIVEE SII PAR RAPPORT A SIGMA =
!          SIJ/SII*(K_IK*K_KL-1/3*K_IJ*K_KL)
    do i = 1, ndt
        dsiids(i) = zero
        do j = 1, ndt
            dsdsig(j, i) = dikdjl(j, i)-dijdkl(j, i)
            dsiids(i) = dsiids(i)+devsig(j)*dsdsig(j, i)/sii
        end do
    end do
!
! --- CALCUL DE DHDS*DSIIDS
    call lcprte(dhds, dsiids, mat1)
!
! --- CALCUL DE D2HDS2
    call lkd2hs(nmat, materf, devsig, sii, rcos3t, &
                dhds, d2hds2)
!
! --- CALCUL DE D2HDSIGMA
    d2hdsi(1:ndt, 1:ndt) = matmul(d2hds2(1:ndt, 1:ndt), dsdsig(1:ndt, 1:ndt))
!
! --- CONSTRUCTION DE SII*D2HDSDSIGMA = MAT2
    mat2(1:ndt, 1:ndt) = sii*d2hdsi(1:ndt, 1:ndt)
!
! --- ADDITION DE MAT2 + MAT1 = MAT3
    mat3(1:ndt, 1:ndt) = mat1(1:ndt, 1:ndt)+mat2(1:ndt, 1:ndt)
!
! --- CALCUL DE COEFH*MAT3 = MAT2
    mat2(1:ndt, 1:ndt) = coefh*mat3(1:ndt, 1:ndt)
!
! --- CONSTRUCTION DE DHTDSIGMA = DHTDS*DSDSIG
    mat1(1:ndt, 1:ndt) = coefh*dsdsig(1:ndt, 1:ndt)
    do i = 1, ndt
        dhtds(i) = zero
        do j = 1, ndt
            dhtds(i) = dhtds(i)+dhds(j)*mat1(j, i)/sii
        end do
    end do
!
! --- CONSTRUCTION PRODUIT TENSORIEL DE MAT1 = DEVSIG*DHTDSIGMA
    call lcprte(devsig, dhtds, mat1)
!
! --- CALCUL DE HTHETA/SII*DSDSIG = MAT3
    mat3(1:ndt, 1:ndt) = (varh(3)/sii)*dsdsig(1:ndt, 1:ndt)
!
! --- MAT5 = HTHETA*DEVSIG*DSIIDS/SII**2
    call lcprte(devsig, dsiids, mat4)
    mat5(1:ndt, 1:ndt) = (varh(3)/sii**2)*mat4(1:ndt, 1:ndt)
!
! --- MAT4 = MAT2+MAT1+MAT3-MAT5
    do i = 1, ndt
        do j = 1, ndt
            mat4(i, j) = mat2(i, j)+mat1(i, j)+mat3(i, j)-mat5(i, j)
        end do
    end do
!
! --- D2SHDS = DSDSIG.MAT4
    d2shds(1:ndt, 1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), mat4(1:ndt, 1:ndt))
!
end subroutine
