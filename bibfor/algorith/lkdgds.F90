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
subroutine lkdgds(nmat, materf, para, vara, devsig, &
                  i1, val, ds2hds, vecn, dfds, &
                  bprimp, nvi, vint, dhds, dgds, &
                  iret)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!     ------------------------------------------------------------------
!     CALCUL DE DERIVEE DE G PAR RAPPORT A SIGMA
!     POUR G PLASTIQUE OU VISQUEUX
!     IN  NMAT   : DIMENSION TABLE PARAMETRES MATERIAU
!         MATERF : TABLE DES PARAMETRES MATERIAU A T+DT
!         PARA   : VECTEUR CONTENANT AXI, SXI ET MXI
!         VARA   : VECTEUR CONTENANT ADXI, BDXI, DDXI ET KDXI
!         DEVSIG : DEVIATEUR DES CONTRAINTES
!         I1     : TRACE DES CONTRAINTES
!         VAL    : BOOLEEN PRECISANT DILATANCE EN PRE(0) OU POST-PIC(1)
!         VECN   : VECTEUR UNITAIRE DE DILATANCE
!         DFDS   : DERIVEE SEUIL PAR RAPPORT A SIGMA
!         DS2HDS : DERIVEE DE SII*H PAR RAPPORT A SIGMA
!         BPRIMP : PARAMETRE DE DILATANCE FCTN SIGMA
!         NVI    : NOMBRE DE VARIABLES INTERNES
!         VINT   : VARIABLES INTERNES
!         DHDS   : DERIVEE DE H(THETA) PAR RAPPORT A SIGMA
!     OUT DGDS   : DERIVEE DU POTENTIEL G PAR RAPPORT A SIGMA (NDTXNDT)
!         IRET   : CODE RETOUR
!     ------------------------------------------------------------------
#include "asterfort/cos3t.h"
#include "asterfort/lcprte.h"
#include "asterfort/lkd2fs.h"
#include "asterfort/lkd2sh.h"
#include "asterfort/lkdnds.h"
#include "asterfort/lkhtet.h"
    integer(kind=8) :: iret, nmat, nvi, val
    real(kind=8) :: materf(nmat, 2), dgds(6, 6), vecn(6), dfds(6), vint(nvi)
    real(kind=8) :: para(3), vara(4), ds2hds(6), devsig(6), i1, bprimp
    real(kind=8) :: dhds(6)
!
    integer(kind=8) :: i, j, ndi, ndt
    real(kind=8) :: d2fds2(6, 6), d2fdsn(6)
    real(kind=8) :: dfdnpn(6, 6), dfpndn(6, 6)
    real(kind=8) :: zero, dfdsdn(6), dfdsvn, d2shds(6, 6), varh(3)
    real(kind=8) :: h0e, h0c, htheta, rcos3t, lgleps, patm
    real(kind=8) :: dndsig(6, 6), d2fn2(6, 6)
    parameter(zero=0.d0)
    parameter(lgleps=1.0d-8)
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!     ------------------------------------------------------------------
!
    patm = materf(1, 2)
!
! --- APPEL A HOC ET  H(THETA)
    rcos3t = cos3t(devsig, patm, lgleps)
    call lkhtet(nmat, materf, rcos3t, h0e, h0c, &
                htheta)
    varh(1) = h0e
    varh(2) = h0c
    varh(3) = htheta
!
! --- CONSTRUCTION D2FDS2
    call lkd2sh(nmat, materf, varh, dhds, devsig, &
                rcos3t, d2shds, iret)
!
    call lkd2fs(nmat, materf, para, vara, varh, &
                i1, devsig, ds2hds, d2shds, d2fds2, &
                iret)
!
! --- CONSTRUCTION DNDSIG
    call lkdnds(nmat, materf, i1, devsig, bprimp, &
                nvi, vint, val, para, dndsig)
!
! --- CONSTRUCTION DE (D2FDS2:N)*N
    do i = 1, ndt
        d2fdsn(i) = zero
        do j = 1, ndt
            d2fdsn(i) = d2fdsn(i)+vecn(j)*d2fds2(j, i)
        end do
    end do
    call lcprte(vecn, d2fdsn, d2fn2)
!
! --- CONSTRUCTION DE (DFDS*DNDSIG)*VECN
    do i = 1, ndt
        dfdsdn(i) = zero
        do j = 1, ndt
            dfdsdn(i) = dfdsdn(i)+dfds(j)*dndsig(j, i)
        end do
    end do
    call lcprte(vecn, dfdsdn, dfdnpn)
!
! --- CONSTRUCTION DE (DFDS:VECN)*DNDSIG
    dfdsvn = dot_product(dfds(1:ndt), vecn(1:ndt))
    dfpndn(1:ndt, 1:ndt) = dfdsvn*dndsig(1:ndt, 1:ndt)
!
! --- CONSTRUCTION DE DGDS
    do i = 1, ndt
        do j = 1, ndt
            dgds(i, j) = d2fds2(i, j)-d2fn2(i, j)-dfdnpn(i, j)-dfpndn(i, j)
        end do
    end do
!
end subroutine
