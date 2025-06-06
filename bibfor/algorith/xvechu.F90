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

subroutine xvechu(ndim, nnop, nnops, ddls, ddlm, pla, &
                  lamb, am, delta, r, p, ffp, jac, ffc, vect, &
                  ncompn, jheavn, ifiss, nfiss, nfh, &
                  ifa, jheafa, ncomph)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/hmdeca.h"
#include "asterfort/transp.h"
#include "asterfort/prmave.h"
#include "asterfort/xcalc_saut.h"
#include "asterfort/xcalc_code.h"
!
! ======================================================================
! person_in_charge: daniele.colombo at ifpen.fr
!
!
! ROUTINE MODELE HM-XFEM
!
! CALCUL DES SECONDS MEMBRES VECT (EQUILIBRE MECANIQUE + LOI INTERFACE)
!
! ----------------------------------------------------------------------
!
    integer :: i, j, k, ndim, ier, nnop, ddls, ddlm, nnops, in, pla(27)
    integer :: pli, ncompn, jheavn, nfiss, hea_fa(2)
    integer :: ifiss, nfh, jheafa, ifa, ncomph, ifh, dec
    real(kind=8) :: h(3), hfix(3), ptr(3, 3), lamb(3), am(3), delta(6)
    real(kind=8) :: r, p(3, 3), ffp(27), jac, ffc(16), ffi
    real(kind=8) :: vect(560), coefi
    aster_logical :: lmultc
!
!   INITIALISATIONS
    lmultc = nfiss .gt. 1
    h(:) = 0.d0
    hfix(:) = 0.d0
    ptr(:, :) = 0.d0
!
    do i = 1, ndim
        h(i) = -lamb(i)-r*am(i)+r*delta(i)
    end do
!
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    else
        hea_fa(1) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+1)
        hea_fa(2) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+2)
    end if
!
!   CONVERSION DE H EN BASE FIXE : {HFIX} = [P]T {H}
!
    call transp(p, 3, ndim, ndim, ptr, 3)
!
    call prmave(0, ptr, 3, ndim, ndim, &
                h, ndim, hfix, ndim, ier)
!
    coefi = xcalc_saut(1, 0, 1)
    do i = 1, nnop
        call hmdeca(i, ddls, ddlm, nnops, in, dec)
        do ifh = 1, nfh
            coefi = xcalc_saut(zi(jheavn-1+ncompn*(i-1)+ifh), &
                               hea_fa(1), &
                               hea_fa(2), &
                               zi(jheavn-1+ncompn*(i-1)+ncompn))
            do j = 1, ndim
                vect(in+(ndim+dec)*ifh+j) = vect(in+(ndim+dec)*ifh+j) &
                                            -coefi*ffp(i)*hfix(j)*jac
            end do
        end do
    end do
!
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
        do k = 1, ndim
            vect(pli+2+k) = vect(pli+2+k)+(am(k)-delta(k))*ffi*jac
        end do
    end do
!
end subroutine
