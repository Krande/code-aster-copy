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
subroutine lkd2fs(nmat, materf, para, vara, varh, &
                  i1, devsig, ds2hds, d2shds, d2fds2, &
                  iret)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!     ------------------------------------------------------------------
!     CALCUL DE DERIVEE 2NDE DE F PAR RAPPORT A SIGMA
!     IN  NMAT   : DIMENSION TABLE DES PARAMETRES MATERIAU
!         MATERF : PARAMETRES MATERIAU A T+DT
!         I1     : TRACE DES CONTRAINTES
!         DEVSIG : DEVIATEUR DES CONTRAINTES
!         PARA   : PARAMETRES AXI, SXI, MXI
!         VARA   : PARAMTERES ADXI,BDXI,DDXI,KDXI
!         VARH   : VECTEUR CONTENANT H0E,H0C ET HTHETA
!         DS2HDS : DERIVEE DE SII*H PAR RAPPORT A SIGMA
!         D2SHDS : DERIVVE 2NDE DE SII*H PAR RAPPORT A SIGMA
!     OUT D2FDS2 :  DERIVEE 2NDE F PAR RAPPORT A SIGMA (NDT X NDT)
!         IRET   :  CODE RETOUR
!     ------------------------------------------------------------------
#include "asterfort/lcprte.h"
    integer(kind=8) :: iret, nmat
    real(kind=8) :: d2fds2(6, 6), para(3), vara(4), materf(nmat, 2)
    real(kind=8) :: devsig(6), i1, ds2hds(6), varh(3), d2shds(6, 6)
!
    integer(kind=8) :: ndi, ndt, i
    real(kind=8) :: sigc, sii, coef1, coef2, vident(6), zero, un, vect1(6)
    real(kind=8) :: mat1(6, 6), mat2(6, 6), mat3(6, 6), ucri, deux
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!     ------------------------------------------------------------------
!
! --- RECUPERATION PARAMETRES MATERIAU
    sigc = materf(3, 2)
!
! --- CONSTRUCTION DE SII
    sii = norm2(devsig(1:ndt))
!
! --- CONSTRUCTION COEF1 = A*SIGC*H0C*(A-1)(AD*SII*H+B*I1+D)^(A-2)
    ucri = vara(1)*sii*varh(3)+vara(2)*i1+vara(3)
    if (ucri .le. zero) then
        ucri = zero
        coef1 = zero
        coef2 = un
    else
        coef1 = para(1)*sigc*varh(2)*(para(1)-un)*ucri**(para(1)-deux)
! --- CONSTRUCTION COEF2 = A*SIGC*H0C(AD*SII*H+B*I1+D)^(A-1)
        coef2 = un-(vara(1)*para(1)*sigc*varh(2)*ucri**(para(1)-un))
    end if
!
! --- CONSTRUCTION VECTEUR IDENTITE
    vident(:) = zero
    do i = 1, ndi
        vident(i) = un
    end do
!
! --- CONSTRUCTION (A*DS2HDS+B*VIDENT)
    do i = 1, ndt
        vect1(i) = vara(1)*ds2hds(i)+vara(2)*vident(i)
    end do
! --- CONSTRUCTION PRODUIT TENSORIEL COEF1*(VECT1 X VECT1)
    call lcprte(vect1, vect1, mat1)
    mat2(1:ndt, 1:ndt) = coef1*mat1(1:ndt, 1:ndt)
!
! --- CONSTRUCTION PRODUIT COEF2*D2SHDS
    mat3(1:ndt, 1:ndt) = coef2*d2shds(1:ndt, 1:ndt)
!
! --- CONSTRUCTION DIFFERENCE MAT3-MAT2 = D2FDS2
    d2fds2(1:ndt, 1:ndt) = mat3(1:ndt, 1:ndt)-mat2(1:ndt, 1:ndt)
!
end subroutine
