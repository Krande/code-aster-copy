! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine xhmco4(ndim, nnop, nnops, pla, nd, tau1, &
                  tau2, ffc, nddls, jac, ffp, &
                  nddlm, mmat, ifiss, nfiss, nfh, &
                  ifa, jheafa, ncomph, jheavn, ncompn)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/hmdeca.h"
#include "asterfort/transp.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nnop, nnops, nddls, pla(27)
    integer :: nddlm
    real(kind=8) :: mmat(560, 560)
    real(kind=8) :: ffp(27), jac, ffc(16)
    real(kind=8) :: nd(3), tau1(3), tau2(3)
    integer :: ifiss, nfiss, nfh, ifa, jheafa, ncomph
    integer :: jheavn, ncompn
! ======================================================================
! person_in_charge: daniele.colombo at ifpen.fr
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES MATRICES DE COHESION, LOI CZM_LIN_MIX
! --- PARTIE INDEPENDANTE DE LA LOI D'INTERFACE
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS SOMMETS DE L'ELEMENT DE REF PARENT
! IN  PLA    : PLACE DES LAGRANGES DANS LA NUMEROTATION
! IN  ND     : DIRECTION NORMALE
! IN  TAU1   : DIRECTION TANGENTE 1
! IN  TAU2   : DIRECTION TANGENTE 2
! IN  NDDLS  : NOMBRE DE DDLS DES NOEUDS SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  NDDLM  : NOMBRE DE DDLS DES NOEUDS MILIEU
! I/O MMAT   : MATRICE ELEMENTAITRE DE COHESION
!
    integer :: i, j, k, l, pli, plj, jn, hea_fa(2), dec, ifh
    real(kind=8) :: dside2(3, 3), ptr(3, 3), temp(3, 3), au(3, 3)
    real(kind=8) :: p(3, 3), coefj, ffi, ffj
    aster_logical :: lmultc
!
! ----------------------------------------------------------------------
!
!     INITIALISATION
    au(:, :) = 0.d0
    dside2(:, :) = 0.d0
    temp(:, :) = 0.d0
    ptr(:, :) = 0.d0
    p(:, :) = 0.d0
    lmultc = nfiss .gt. 1
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    else
        hea_fa(1) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+1)
        hea_fa(2) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+2)
    end if
!
! idem, il va falloir introduire les matrices de passage
    do i = 1, ndim
        p(1, i) = nd(i)
    end do
    do i = 1, ndim
        p(2, i) = tau1(i)
    end do
    if (ndim .eq. 3) then
        do i = 1, ndim
            p(3, i) = tau2(i)
        end do
    end if
!
! on construit la transposee de la matrice de passage
    call transp(p, 3, ndim, ndim, ptr, &
                3)
!
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
!
        do j = 1, nnop
            call hmdeca(j, nddls, nddlm, nnops, jn, dec)
!
            do ifh = 1, nfh
                coefj = xcalc_saut(zi(jheavn-1+ncompn*(j-1)+ifh), &
                                   hea_fa(1), &
                                   hea_fa(2), &
                                   zi(jheavn-1+ncompn*(j-1)+ncompn))
                do l = 1, ndim
                    do k = 1, ndim
! on remplit A : matrice [u*] / mu
                        mmat(pli-1+3+2*ndim+k, jn+(ndim+dec)*ifh+l) = &
                            mmat(pli-1+3+2*ndim+k, jn+(ndim+dec)*ifh+l)+ &
                            coefj*ffi*p(k, l)*ffp(j)*jac
! et sa transposee
                        mmat(jn+(ndim+dec)*ifh+l, pli-1+3+2*ndim+k) = &
                            mmat(jn+(ndim+dec)*ifh+l, pli-1+3+2*ndim+k)+ &
                            coefj*ffi*p(k, l)*ffp(j)*jac
                    end do
                end do
            end do
        end do
    end do
!
! on remplit B : matrice w* / mu
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
!
        do j = 1, nnops
            plj = pla(j)
            ffj = ffc(j)
            do l = 1, ndim
! on remplit B
                mmat(pli-1+3+ndim+l, plj-1+3+2*ndim+l) = mmat(pli-1+3+ndim+l, plj-1+3+2*ndim+l)- &
                                                         ffi*ffj*jac
! et sa transposee
                mmat(plj-1+3+2*ndim+l, pli-1+3+ndim+l) = mmat(plj-1+3+2*ndim+l, pli-1+3+ndim+l)- &
                                                         ffi*ffj*jac
            end do
! on remplit  w* / pf
            mmat(pli-1+3+ndim+1, plj) = mmat(pli-1+3+ndim+1, plj)-ffi*ffj*jac
        end do
    end do
!
end subroutine
