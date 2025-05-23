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

subroutine xmmatu(ndim, nnop, nnops, ddls, ddlm, pla, &
                  dsidep, p, r, ffp, jac, ffc, nd, mmat, &
                  jheavn, ncompn, ifiss, nfiss, nfh, ifa, &
                  jheafa, ncomph)

    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/hmdeca.h"
#include "asterfort/promat.h"
#include "asterfort/transp.h"
#include "asterfort/prmave.h"
#include "asterfort/xcalc_saut.h"
#include "asterfort/xcalc_code.h"
!
! person_in_charge: daniele.colombo at ifpen.fr
! ======================================================================
!
! ROUTINE MODELE HM-XFEM (CAS DE LA FRACTURE)
!
! CALCUL DE LA MATRICE MMAT
!
! ----------------------------------------------------------------------
!
    integer :: i, j, ndim, ier1, ier2, nnop, ddls, ddlm, nnops, in, jn
    integer :: k, l, pla(27), pli, plj, jheavn, ncompn, nfiss, hea_fa(2)
    integer :: ifiss, nfh, ifa, ncomph, jheafa, ifh, jfh, dec, dej
    real(kind=8) :: unity(3, 3), dside2(3, 3), alocal(3, 3), ptr(3, 3), pdotal(3, 3)
    real(kind=8) :: kdotal(3, 3), au(3, 3), dsidep(6, 6), temp(3, 3), knd(3), knloc(3)
    real(kind=8) :: r, p(3, 3), ffp(27), jac, ffj, ffc(16), ffi, nd(3)
    real(kind=8) :: mmat(560, 560), dside3(3, 3), coefi, coefj
    aster_logical :: lmultc
!
!   INITIALISATION
    lmultc = nfiss .gt. 1
    unity(:, :) = 0.d0
    alocal(:, :) = 0.d0
    ptr(:, :) = 0.d0
    pdotal(:, :) = 0.d0
    kdotal(:, :) = 0.d0
    au(:, :) = 0.d0
    dside2(:, :) = 0.d0
    dside3(:, :) = 0.d0
    temp(:, :) = 0.d0
    knd(:) = 0.d0
    knloc(:) = 0.d0
!
!   MATRICE -ID+R DSIDEP
!
    do i = 1, ndim
        unity(i, i) = 1.d0
    end do
!
    do i = 1, ndim
        do j = 1, ndim
            dside2(i, j) = dsidep(i, j)
            dside3(i, j) = -dsidep(i, j)
            alocal(i, j) = -unity(i, j)+r*dside2(i, j)
        end do
    end do
!
!   MATRICE [P]T[ALOCAL]
!
    call transp(p, 3, ndim, ndim, ptr, 3)
!
    call promat(ptr, 3, ndim, ndim, alocal, &
                3, ndim, ndim, pdotal)
!
!   MATRICE [P]T[DSIDEP]
!
    call promat(ptr, 3, ndim, ndim, dside3, &
                3, ndim, ndim, temp)
!
    call promat(temp, 3, ndim, ndim, p, &
                3, ndim, ndim, kdotal)
!
    call prmave(0, kdotal, 3, ndim, ndim, &
                nd, ndim, knd, ndim, ier1)
!
!   VECTEUR [DSIDEP].ND
!
    temp(:, :) = 0.d0
    call promat(dside2, 3, ndim, ndim, p, &
                3, ndim, ndim, temp)
!
    call prmave(0, temp, 3, ndim, ndim, &
                nd, ndim, knloc, ndim, ier2)
!
!   MATRICE TANGENTE EN BASE FIXE [P]T [DSIDEP] [P]
!
    temp(:, :) = 0.d0
    call promat(ptr, 3, ndim, ndim, alocal, &
                3, ndim, ndim, temp)
!
    call promat(temp, 3, ndim, ndim, p, &
                3, ndim, ndim, au)
!
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    else
        hea_fa(1) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+1)
        hea_fa(2) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+2)
    end if
!
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
        do k = 1, ndim
            do j = 1, nnop
                call hmdeca(j, ddls, ddlm, nnops, jn, dec)
                do ifh = 1, nfh
                    coefj = xcalc_saut(zi(jheavn-1+ncompn*(j-1)+ifh), &
                                       hea_fa(1), &
                                       hea_fa(2), &
                                       zi(jheavn-1+ncompn*(j-1)+ncompn))
                    do l = 1, ndim
! INDICES INVERSES MATRICE INTERFACE
                        mmat(pli+2+k, jn+(ndim+dec)*ifh+l) = &
                            mmat(pli+2+k, jn+(ndim+dec)*ifh+l)- &
                            coefj*ffp(j)*pdotal(l, k)*ffi*jac
! INDICES MEME ORDRE MATRICE EQUILIBRE
                        mmat(jn+(ndim+dec)*ifh+l, pli+2+k) = &
                            mmat(jn+(ndim+dec)*ifh+l, pli+2+k)- &
                            coefj*ffp(j)*pdotal(l, k)*ffi*jac
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, nnop
        call hmdeca(i, ddls, ddlm, nnops, in, dec)
        do ifh = 1, nfh
            coefi = xcalc_saut(zi(jheavn-1+ncompn*(i-1)+ifh), &
                               hea_fa(1), &
                               hea_fa(2), &
                               zi(jheavn-1+ncompn*(i-1)+ncompn))
            do j = 1, nnop
                call hmdeca(j, ddls, ddlm, nnops, jn, dej)
                do jfh = 1, nfh
                    coefj = xcalc_saut(zi(jheavn-1+ncompn*(j-1)+jfh), &
                                       hea_fa(1), &
                                       hea_fa(2), &
                                       zi(jheavn-1+ncompn*(j-1)+ncompn))
                    do k = 1, ndim
                        do l = 1, ndim
                            mmat(in+(ndim+dec)*ifh+k, jn+(ndim+dej)*jfh+l) = &
                                mmat(in+(ndim+dec)*ifh+k, jn+(ndim+dej)*jfh+l)- &
                                r*au(k, l)*coefi*ffp(i)*coefj*ffp(j)*jac
                        end do
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, nnop
        call hmdeca(i, ddls, ddlm, nnops, in, dec)
        do ifh = 1, nfh
            coefi = xcalc_saut(zi(jheavn-1+ncompn*(i-1)+ifh), &
                               hea_fa(1), &
                               hea_fa(2), &
                               zi(jheavn-1+ncompn*(i-1)+ncompn))
            do k = 1, ndim
                do j = 1, nnops
                    plj = pla(j)
                    ffj = ffc(j)
                    mmat(in+(ndim+dec)*ifh+k, plj) = &
                        mmat(in+(ndim+dec)*ifh+k, plj)+coefi*r*ffp(i)*knd(k)*ffj*jac
                end do
            end do
        end do
    end do
!
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
        do k = 1, ndim
            do j = 1, nnops
                plj = pla(j)
                ffj = ffc(j)
                do l = 1, ndim
                    mmat(pli+2+k, plj+2+l) = &
                        mmat(pli+2+k, plj+2+l)-ffi*dside2(k, l)*ffj*jac
                end do
            end do
        end do
    end do
!
    do i = 1, nnops
        pli = pla(i)
        ffi = ffc(i)
        do k = 1, ndim
            do j = 1, nnops
                plj = pla(j)
                ffj = ffc(j)
                mmat(pli+2+k, plj) = mmat(pli+2+k, plj)- &
                                     ffi*knloc(k)*ffj*jac
            end do
        end do
    end do
!
end subroutine
