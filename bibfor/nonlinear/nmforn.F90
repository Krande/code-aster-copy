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
subroutine nmforn(ndim, nno1, nno2, npg, iw, &
                  vff1, vff2, idfde1, geom, vect)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmgvdn.h"
#include "asterfort/nmmabu.h"
#include "asterfort/r8inir.h"
#include "asterfort/terefe.h"
    integer(kind=8) :: ndim, nno1, nno2, npg, idfde1, iw
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    real(kind=8) :: geom(ndim, nno1)
    real(kind=8) :: vect(*)
! ---------------------------------------------------------------------
!
!     FORC_NODA POUR GRAD_VARI (2D ET 3D)
!
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  IDFDE1  : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE E)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE E)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS DE REFERENCE (INDICE)
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! ---------------------------------------------------------------------
!
    aster_logical :: grand, axi
    integer(kind=8) :: nddl, ndimsi, g, n, i, kl, kk
    integer(kind=8) :: iu(3*27), ia(8)
    real(kind=8) :: dfdi1(27, 3)
    real(kind=8) :: r, wg, b(6, 3, 27)
    real(kind=8) :: t1, sigref, varref
! ---------------------------------------------------------------------
!
! - INITIALISATION
!
    grand = .false.
    axi = .false.
    nddl = nno1*ndim+nno2
    ndimsi = 2*ndim
!
!
! --- VALEURS DE REFERENCE POUR REFE_FORC_NODA
!
    call terefe('SIGM_REFE', 'MECA_GRADVARI', sigref)
    call terefe('VARI_REFE', 'MECA_GRADVARI', varref)
!
    call r8inir(nddl, 0.d0, vect, 1)
!
    call nmgvdn(ndim, nno1, nno2, iu, ia)
!
    do g = 1, npg
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF POUR U
!
        call dfdmip(ndim, nno1, axi, geom, g, &
                    iw, vff1(1, g), idfde1, r, wg, &
                    dfdi1)
        call nmmabu(ndim, nno1, axi, grand, dfdi1, &
                    b)
        do n = 1, nno1
            do i = 1, ndim
                kk = iu(nno1*(i-1)+n)
                t1 = 0
                do kl = 1, ndimsi
                    t1 = t1+abs(b(kl, i, n))
                end do
                vect(kk) = vect(kk)+wg*t1*sigref
            end do
        end do
!
        do n = 1, nno2
            kk = ia(n)
            vect(kk) = vect(kk)+wg*vff2(n, g)*varref
        end do
!
!
    end do
!
end subroutine
