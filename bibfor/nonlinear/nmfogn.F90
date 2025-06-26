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
subroutine nmfogn(ndim, nno1, nno2, npg, iw, &
                  vff1, vff2, idfde1, idfde2, geom, &
                  typmod, mat, ddl, sigm, vect)
!
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmgvdn.h"
#include "asterfort/nmmabu.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
    character(len=8) :: typmod(*)
    integer(kind=8) :: ndim, nno1, nno2, npg, idfde1, idfde2, iw, mat
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    real(kind=8) :: geom(ndim, nno1)
    real(kind=8) :: sigm(2*ndim+1, npg), vect(*), ddl(*)
! ---------------------------------------------------------------------
!
!     FORC_NODA POUR GVNO (2D ET 3D)
!
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  IDFDE1  : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE E)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE E)
! IN  IDFDE2  : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE E)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS DE REFERENCE (INDICE)
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODEELISATION
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! MEM DFDI2   :
! ---------------------------------------------------------------------
!
    integer(kind=8) :: k2(1)
    character(len=8) :: fami, poum
    character(len=16) :: nom(1)
!
    aster_logical :: grand, axi, nax
    integer(kind=8) :: nddl, ndimsi, g, n, i, kl, kk
    integer(kind=8) :: iu(3*27), ia(8), kpg, spt
    real(kind=8) :: rac2, c, val(1)
    real(kind=8) :: dfdi1(27, 3)
    real(kind=8) :: av, ag(3), bp
    real(kind=8) :: r, wg, b(6, 3, 27)
    real(kind=8) :: sigma(6), t1, t2
    real(kind=8) :: dfdi2(8*3)
! ---------------------------------------------------------------------
!
!
!
! - INITIALISATION
!
    rac2 = sqrt(2.d0)
    grand = .false.
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2
    ndimsi = 2*ndim
!
    call r8inir(nddl, 0.d0, vect, 1)
!
    nom(1) = 'C_GRAD_VARI'
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(fami, kpg, spt, poum, mat, &
                ' ', 'NON_LOCAL', 0, ' ', [0.d0], &
                1, nom, val, k2, 2)
    c = val(1)
!
    call nmgvdn(ndim, nno1, nno2, iu, ia)
!
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
!
    do g = 1, npg
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF POUR U
!
        call dfdmip(ndim, nno1, axi, geom, g, &
                    iw, vff1(1, g), idfde1, r, wg, &
                    dfdi1)
        nax = .false.
        call nmmabu(ndim, nno1, nax, grand, dfdi1, &
                    b)
        if (axi) then
            do n = 1, nno1
                b(3, 1, n) = vff1(n, g)/r
            end do
        end if
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF POUR A
!
        call dfdmip(ndim, nno2, axi, geom, g, &
                    iw, vff2(1, g), idfde2, r, wg, &
                    dfdi2)
!
        av = 0
        do n = 1, nno2
            av = av+vff2(n, g)*ddl(ia(n))
        end do
!
        do i = 1, ndim
            ag(i) = 0
            do n = 1, nno2
                ag(i) = ag(i)+dfdi2(nno2*(i-1)+n)*ddl(ia(n))
            end do
        end do
!
        do kl = 1, 3
            sigma(kl) = sigm(kl, g)
        end do
        do kl = 4, ndimsi
            sigma(kl) = sigm(kl, g)*rac2
        end do
        bp = sigm(ndimsi+1, g)
!
!      VECTEUR FINT:U
!
        do n = 1, nno1
            do i = 1, ndim
                kk = iu(nno1*(i-1)+n)
                t1 = 0
                do kl = 1, ndimsi
                    t1 = t1+sigma(kl)*b(kl, i, n)
                end do
                vect(kk) = vect(kk)+wg*t1
            end do
        end do
!
!      VECTEUR FINT:A
!
        do n = 1, nno2
            t1 = vff2(n, g)*bp
            t2 = 0
            do i = 1, ndim
                t2 = t2+c*dfdi2(nno2*(i-1)+n)*ag(i)
            end do
            kk = ia(n)
            vect(kk) = vect(kk)+wg*(t2+t1)
        end do
!
    end do
!
!
end subroutine
