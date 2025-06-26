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
subroutine lkdnds(nmat, materf, i1, devsig, bprimp, &
                  nvi, vint, val, para, dndsig)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!     ------------------------------------------------------------------
!     CALCUL DE DERIVEE DE N PAR RAPPORT A SIGMA
!     IN  NMAT   : DIMENSION TABLE DES PARAMETRES MATERIAU
!         MATERF : PARAMETRES MATERIAU A T+DT
!         DEVSIG : DEVIATEUR DES CONTRAINTES
!         I1     : TRACE DES CONTRAINTES
!         BPRIMP : PARAMETRE DE DILATANCE FCTN SIGMA
!         NVI    : NOMBRE DE VARIABLES INTERNES
!         VINT   : VARIABLES INTERNES
!         VAL    : BOOLEEN PRECISANT DILATANCE EN PRE(0) OU POST-PIC(1)
!         PARA   : VECTEUR CONTENANT AXI, SXI ET MXI
!
!     OUT DNDISG :  DERIVEE DE N PAR RAPPORT A SIGMA (NDT X NDT)
!     ------------------------------------------------------------------
#include "asterfort/lcprte.h"
#include "asterfort/lkdbds.h"
    integer(kind=8) :: nmat, nvi, val
    real(kind=8) :: materf(nmat, 2), dndsig(6, 6), devsig(6), i1
    real(kind=8) :: bprimp, vint(nvi), para(3)
!
    integer(kind=8) :: ndt, ndi, i, j
    real(kind=8) :: zero, un
    real(kind=8) :: dsdsig(6, 6), di1dsi(6), sii
    real(kind=8) :: dbetds(6), dbetdi, mident(6, 6)
    real(kind=8) :: dsiids(6, 6), kron2(6, 6)
    real(kind=8) :: unstro, unssii, trois, kron(6), dbdsig(6)
    real(kind=8) :: devbds(6, 6), dsidsi(6, 6), sdsids(6, 6), didbds(6, 6)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(trois=3.d0)
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!     ------------------------------------------------------------------
!
! ----------------------------
! --- 1) CALCUL TERMES COMMUNS
! ----------------------------
! --- CONSTRUCTION DE SII
    sii = norm2(devsig(1:ndt))
!
! --- CONSTRUCTION DE KRONECKER
    kron(:) = zero
    do i = 1, ndi
        kron(i) = un
    end do
!
! --- CONSTRUCTION DE MATRICE IDENTITE
    mident(:, :) = zero
    do i = 1, ndt
        mident(i, i) = un
    end do
!
! --- CONSTRUCTION DE DI1/DSIGMA
    di1dsi(:) = zero
    do i = 1, ndi
        di1dsi(i) = un
    end do
!
! --- CONSTRUCTION DE DS/DSIGMA
    unstro = un/trois
    call lcprte(kron, kron, kron2)
    kron2(1:ndt, 1:ndt) = unstro*kron2(1:ndt, 1:ndt)
    dsdsig(1:ndt, 1:ndt) = mident(1:ndt, 1:ndt)-kron2(1:ndt, 1:ndt)
!
! --- CONSTRUCTION DE DSII/D(DEVSIG)
    unssii = un/sii
    dsiids(1:ndt, 1) = unssii*devsig(1:ndt)
!
! --- CALCUL DE D(BPRIME)/D(DEVSIG) ET D(BPRIME)/D(I1)
    call lkdbds(nmat, materf, i1, devsig, nvi, &
                vint, para, val, dbetds, dbetdi)
!
! --- CONSTRUCTION DE D(BPRIME)/DSIGMA
    do i = 1, ndt
        dbdsig(i) = zero
        do j = 1, ndt
            dbdsig(i) = dbdsig(i)+dbetds(j)*dsdsig(j, i)
        end do
        dbdsig(i) = dbdsig(i)+dbetdi*di1dsi(i)
    end do
!
! --- PRODUIT TENSORIEL DE DEVSIG*D(BPRIME)/DSIGMA
    call lcprte(devsig, dbdsig, devbds)
!
! --- CALCUL DE DEVSIG*DSII/DSIGMA
    dsidsi(1:ndt, 1) = matmul(dsdsig(1:ndt, 1:ndt), dsiids(1:ndt, 1))
    call lcprte(devsig, dsidsi, sdsids)
!
! --- PRODUIT TENSORIEL DE KRON*DBPRIME/DSIGMA
    call lcprte(kron, dbdsig, didbds)
!
! ------------------------------
! --- 2) ASSEMBLAGE DE DN/DSIGMA
! ------------------------------
    do i = 1, ndt
        do j = 1, ndt
            dndsig(i, j) = ( &
                          ( &
                          devbds(i, j)/sii+bprimp/sii*dsdsig(i, j)-bprimp/sii**2*sdsids(i, j))*sqr&
                          &t(bprimp**2+trois)-bprimp/sqrt(bprimp**2+trois)*(bprimp/sii*devbds(i&
                          &, j)-didbds(i, j)) &
                          )/(bprimp**2+trois &
                          )
        end do
    end do
!
end subroutine
