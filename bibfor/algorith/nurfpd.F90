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
subroutine nurfpd(ndim, nno1, nno2, npg, iw, &
                  vff1, vff2, idff1, vu, vp, &
                  typmod, geomi, sigref, epsref, vect)
! person_in_charge: sebastien.fayolle at edf.fr
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/r8inir.h"
#include "blas/ddot.h"
    integer(kind=8) :: ndim, nno1, nno2, npg, iw, idff1
    integer(kind=8) :: vu(3, 27), vp(27)
    real(kind=8) :: geomi(ndim, nno1)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    real(kind=8) :: sigref, epsref
    real(kind=8) :: vect(*)
    character(len=8) :: typmod(*)
!
!-----------------------------------------------------------------------
!          CALCUL DE REFE_FORC_NODA POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!          ROUTINE APPELEE PAR TE0598
!-----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF1    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  GEOMI   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  SIGREF  : CONTRAINTE DE REFERENCE
! IN  EPSREF  : DEFORMATION DE REFERENCE
! OUT VECT    : REFE_FORC_NODA
!-----------------------------------------------------------------------
!
    aster_logical :: axi
    integer(kind=8) :: nddl, g
    integer(kind=8) :: kl, sa, na, ia, kk
    integer(kind=8) :: ndimsi
    real(kind=8) :: r, w, sigma(6)
    real(kind=8) :: rac2
    real(kind=8) :: f(3, 3)
    real(kind=8) :: def(2*ndim, nno1, ndim)
    real(kind=8) :: t1, dff1(nno1, 4)
    blas_int :: b_incx, b_incy, b_n
!
    data f/1.d0, 0.d0, 0.d0,&
     &                    0.d0, 1.d0, 0.d0,&
     &                    0.d0, 0.d0, 1.d0/
!-----------------------------------------------------------------------
!
! - INITIALISATION
!
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2
    rac2 = sqrt(2.d0)
    ndimsi = 2*ndim
!
    call r8inir(nddl, 0.d0, vect, 1)
!
    do g = 1, npg
!
        call dfdmip(ndim, nno1, axi, geomi, g, &
                    iw, vff1(1, g), idff1, r, w, &
                    dff1)
!
! - CALCUL DE LA MATRICE B EPS_ij=B_ijkl U_kl
! - DEF (XX,YY,ZZ,2/RAC(2)XY,2/RAC(2)XZ,2/RAC(2)YZ)
        if (ndim .eq. 2) then
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = f(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = f(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = 0.d0
                    def(4, na, ia) = (f(ia, 1)*dff1(na, 2)+f(ia, 2)*dff1(na, 1))/rac2
                end do
            end do
!
! - TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                do na = 1, nno1
                    def(3, na, 1) = f(3, 3)*vff1(na, g)/r
                end do
            end if
        else
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = f(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = f(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = f(ia, 3)*dff1(na, 3)
                    def(4, na, ia) = (f(ia, 1)*dff1(na, 2)+f(ia, 2)*dff1(na, 1))/rac2
                    def(5, na, ia) = (f(ia, 1)*dff1(na, 3)+f(ia, 3)*dff1(na, 1))/rac2
                    def(6, na, ia) = (f(ia, 2)*dff1(na, 3)+f(ia, 3)*dff1(na, 2))/rac2
                end do
            end do
        end if
!
! - VECTEUR FINT:U
        do kl = 1, ndimsi
            call r8inir(6, 0.d0, sigma, 1)
            sigma(kl) = sigref
            do na = 1, nno1
                do ia = 1, ndim
                    kk = vu(ia, na)
                    b_n = to_blas_int(2*ndim)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    t1 = ddot(b_n, sigma, b_incx, def(1, na, ia), b_incy)
                    vect(kk) = vect(kk)+abs(w*t1)/ndimsi
                end do
            end do
        end do
!
! - VECTEUR FINT:P
        do sa = 1, nno2
            kk = vp(sa)
            t1 = vff2(sa, g)*epsref
            vect(kk) = vect(kk)+abs(w*t1)
        end do
    end do
end subroutine
