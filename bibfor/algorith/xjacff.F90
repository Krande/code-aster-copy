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
subroutine xjacff(elrefp, elrefc, elc, ndim, fpg, &
                  jinter, ifa, cface, ipg, nnop, &
                  nnops, igeom, jbasec, xg, jac, &
                  ffp, ffpc, dfdi, nd, tau1, &
                  tau2, dfdic)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2b.h"
#include "asterfort/elelin.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "blas/ddot.h"
    integer(kind=8) :: jinter, ifa, cface(30, 6), ipg, nnop, igeom, jbasec, ndim, nnops
    real(kind=8) :: jac, ffp(27), ffpc(27), dfdi(nnop, 3)
    real(kind=8) :: nd(3), tau1(ndim), tau2(ndim), xg(3)
    character(len=8) :: elrefp, fpg, elrefc, elc
    real(kind=8), intent(out), optional :: dfdic(nnops, 3)
!
!
!                   CALCUL DU JACOBIEN DE LA TRANSFORMATION FACETTE
!                       RÉELLE EN 3D À FACETTE DE RÉFÉRENCE 2D
!                   ET DES FF DE L'ÉLÉMENT PARENT AU POINT DE GAUSS
!               ET DE LA NORMALE À LA FACETTE ORIENTÉE DE ESCL -> MAIT
!     ENTREE
!       ELREFP  : TYPE DE L'ELEMENT DE REF PARENT
!       FPG     : FAMILLE DE POINTS DE GAUSS (SCHEMA D'INTEGRATION)
!       PINTER  : COORDONNÉES DES POINTS D'INTERSECTION
!       IFA     : INDINCE DE LA FACETTE COURANTE
!       CFACE   : CONNECTIVITÉ DES NOEUDS DES FACETTES
!       IPG     : NUMÉRO DU POINTS DE GAUSS
!       NNO     : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
!       IGEOM   : COORDONNEES DES NOEUDS DE L'ELEMENT DE REF PARENT
!       IFISS   : FISSURE COURANTE
!
!     SORTIE
!       G       : COORDONNÉES RÉELLES 3D DU POINT DE GAUSS
!       JAC     : PRODUIT DU JACOBIEN ET DU POIDS
!       FF      : FF DE L'ÉLÉMENT PARENT AU POINT DE GAUSS
!       ND      : NORMALE À LA FACETTE ORIENTÉE DE ESCL -> MAIT
!
!     ------------------------------------------------------------------
!
    real(kind=8) :: norme
    real(kind=8) :: grlt(3), grln(3), norm2, ps
    integer(kind=8) :: ibid, nnoc
    integer(kind=8) :: j, k, i, nno, ipoidf, ivff, idfdef, ndimf
    real(kind=8) :: xe(3), coor3d(3*6)
    character(len=8) :: k8bid
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
!
    call elrefe_info(elrefe=elc, fami=fpg, ndim=ndimf, nno=nno, jpoids=ipoidf, &
                     jvf=ivff, jdfde=idfdef)
!
    ASSERT(nno .eq. 3 .or. nno .eq. 6)
    ASSERT(ndim .eq. 3)
!
! --- INITIALISATION
    nd(:) = 0.d0
    grln(:) = 0.d0
    grlt(:) = 0.d0
    tau1(:) = 0.d0
    tau2(:) = 0.d0
    xe(:) = 0.d0
    ffp(:) = 0.d0
    ffpc(:) = 0.d0
    do i = 1, nnop
        do j = 1, ndim
            dfdi(i, j) = 0.d0
        end do
    end do
!
! --- COORDONNÉES DES NOEUDS DE LA FACETTE DANS LE REPERE GLOBAL NDIM
    do i = 1, 18
        coor3d(i) = 0.d0
    end do
    do i = 1, nno
        do j = 1, ndim
            coor3d((i-1)*ndim+j) = zr(jinter-1+ndim*(cface(ifa, i)-1)+j)
        end do
    end do
!     CALCUL DE JAC EN 3D
    k = 2*(ipg-1)*nno
    call dfdm2b(nno, zr(ipoidf-1+ipg), zr(idfdef+k), coor3d, jac, &
                nd)
!
! --- COORDONNEES REELLES 3D DU POINT DE GAUSS IPG
    xg(:) = 0.d0
    do j = 1, nno
        do i = 1, ndim
            xg(i) = xg(i)+zr(ivff-1+nno*(ipg-1)+j)*coor3d(ndim*(j-1)+i)
        end do
    end do
!
! --- CONSTRUCTION DE LA BASE AU POINT DE GAUSS
!
    do j = 1, ndim
        do k = 1, nno
            grln(j) = grln(j)+zr(ivff-1+nno*(ipg-1)+k)*zr(jbasec-1+ndim*ndim*(k-1)+j)
            grlt(j) = grlt(j)+zr(ivff-1+nno*(ipg-1)+k)*zr(jbasec-1+ndim*ndim*(k-1)+j+ndim)
        end do
    end do
!
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ps = ddot(b_n, grln, b_incx, nd, b_incy)
    if (ps .lt. 0.d0) nd(1:3) = -nd(1:3)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ps = ddot(b_n, grlt, b_incx, nd, b_incy)
    do j = 1, ndim
        tau1(j) = grlt(j)-ps*nd(j)
    end do
!
    call normev(tau1, norme)
!
    if (norme .lt. 1.d-12) then
!       ESSAI AVEC LE PROJETE DE OX
        tau1(1) = 1.d0-nd(1)*nd(1)
        tau1(2) = 0.d0-nd(1)*nd(2)
        if (ndim .eq. 3) tau1(3) = 0.d0-nd(1)*nd(3)
        call normev(tau1, norm2)
        if (norm2 .lt. 1.d-12) then
!         ESSAI AVEC LE PROJETE DE OY
            tau1(1) = 0.d0-nd(2)*nd(1)
            tau1(2) = 1.d0-nd(2)*nd(2)
            if (ndim .eq. 3) tau1(3) = 0.d0-nd(2)*nd(3)
            call normev(tau1, norm2)
        end if
        ASSERT(norm2 .gt. 1.d-12)
    end if
    call provec(nd, tau1, tau2)
    call reeref(elrefp, nnop, zr(igeom), xg, 3, &
                xe, ffp, dfdi=dfdi)
!
!
    if (elrefc .eq. elrefp) goto 999
    if (elrefc(1:3) .eq. 'NON') goto 999
!
!     CALCUL DES FF DE L'ÉLÉMENT DE CONTACT EN CE POINT DE GAUSS
    call elelin(3, elrefc, k8bid, nnoc, ibid)
!
    call reeref(elrefc, nnoc, zr(igeom), xg, ndim, &
                xe, ffpc, dfdi=dfdic)
!
999 continue
!
end subroutine
