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
subroutine nifnpd(ndim, nno1, nno2, nno3, npg, &
                  iw, vff1, vff2, vff3, idff1, &
                  vu, vg, vp, typmod, geomi, &
                  sig, ddl, vect)
! person_in_charge: sebastien.fayolle at edf.fr
! aslint: disable=W1306
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmepsi.h"
#include "asterfort/r8inir.h"
#include "blas/ddot.h"
    integer(kind=8) :: ndim, nno1, nno2, nno3, npg, iw, idff1
    integer(kind=8) :: vu(3, 27), vg(27), vp(27)
    real(kind=8) :: geomi(ndim, nno1)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), vff3(nno3, npg)
    real(kind=8) :: sig(2*ndim+1, npg), ddl(*), vect(*)
    character(len=8) :: typmod(*)
!-----------------------------------------------------------------------
!          CALCUL DES FORCES NODALES POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!          ROUTINE APPELEE PAR TE0591
!-----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AU GONFLEMENT
! IN  NNO3    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF1    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES AU GONFLEMENT
! IN  VFF3    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VG      : TABLEAU DES INDICES DES DDL DE GONFLEMENT
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  GEOMI   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  MATE    : MATERIAU CODE
! IN  DDL     : DEGRES DE LIBERTE A L'INSTANT PRECEDENT
! IN  SIG     : CONTRAINTES A L'INSTANT PRECEDENT
! OUT VECT    : FORCES INTERNES
!-----------------------------------------------------------------------
!
    aster_logical :: axi, grand
    integer(kind=8) :: nddl, g
    integer(kind=8) :: sa, ra, na, ia, kk
    real(kind=8) :: deplm(3*27), gonfm(27), gm, r
!    real(kind=8) :: presm(27), pm
    real(kind=8) :: dff1(nno1, ndim)
    real(kind=8) :: fm(3, 3)
    real(kind=8) :: w
    real(kind=8) :: rac2, def(6, 27, 3)
    real(kind=8) :: epsm(6), sigma(6)
    real(kind=8) :: divum
    real(kind=8) :: t1, t2
    blas_int :: b_incx, b_incy, b_n
!
    parameter(grand=.false._1)
!-----------------------------------------------------------------------
!
! - INITIALISATION
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2+nno3
    rac2 = sqrt(2.d0)
!
    call r8inir(nddl, 0.d0, vect, 1)
!
! - EXTRACTION DES CHAMPS
    do na = 1, nno1
        do ia = 1, ndim
            deplm(ia+ndim*(na-1)) = ddl(vu(ia, na))
        end do
    end do
!
    do ra = 1, nno2
        gonfm(ra) = ddl(vg(ra))
    end do
!
!    do sa = 1, nno3
!        presm(sa) = ddl(vp(sa))
!    end do
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
    do g = 1, npg
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
        call r8inir(6, 0.d0, epsm, 1)
        call dfdmip(ndim, nno1, axi, geomi, g, &
                    iw, vff1(1, g), idff1, r, w, &
                    dff1)
        call nmepsi(ndim, nno1, axi, grand, vff1(1, g), &
                    r, dff1, deplm, fm, epsm)
!
        divum = epsm(1)+epsm(2)+epsm(3)
!
! - CALCUL DE LA MATRICE B EPS_ij=B_ijkl U_kl
! - DEF (XX,YY,ZZ,2/RAC(2)XY,2/RAC(2)XZ,2/RAC(2)YZ)
        if (ndim .eq. 2) then
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = 0.d0
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                end do
            end do
!
! - TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                do na = 1, nno1
                    def(3, na, 1) = fm(3, 3)*vff1(na, g)/r
                end do
            end if
        else
            do na = 1, nno1
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = fm(ia, 3)*dff1(na, 3)
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                    def(5, na, ia) = (fm(ia, 1)*dff1(na, 3)+fm(ia, 3)*dff1(na, 1))/rac2
                    def(6, na, ia) = (fm(ia, 2)*dff1(na, 3)+fm(ia, 3)*dff1(na, 2))/rac2
                end do
            end do
        end if
!
! - CALCUL DE LA PRESSION ET DU GONFLEMENT
        b_n = to_blas_int(nno2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        gm = ddot(b_n, vff2(1, g), b_incx, gonfm, b_incy)
!        pm = ddot(nno3,vff3(1,g),1,presm,1)
!
! - CALCUL DES CONTRAINTES MECANIQUES A L'EQUILIBRE
        do ia = 1, 3
            sigma(ia) = sig(ia, g)
        end do
        do ia = 4, 2*ndim
            sigma(ia) = sig(ia, g)*rac2
        end do
!
! - VECTEUR FINT:U
        do na = 1, nno1
            do ia = 1, ndim
                kk = vu(ia, na)
                b_n = to_blas_int(2*ndim)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                t1 = ddot(b_n, sigma, b_incx, def(1, na, ia), b_incy)
                vect(kk) = vect(kk)+w*t1
            end do
        end do
!
! - VECTEUR FINT:G
        t2 = sig(2*ndim+1, g)
        do ra = 1, nno2
            kk = vg(ra)
            t1 = vff2(ra, g)*t2
            vect(kk) = vect(kk)+w*t1
        end do
!
! - VECTEUR FINT:P
        t2 = (divum-gm)
        do sa = 1, nno3
            kk = vp(sa)
            t1 = vff3(sa, g)*t2
            vect(kk) = vect(kk)+w*t1
        end do
    end do
!
end subroutine
