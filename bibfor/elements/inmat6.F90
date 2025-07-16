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
subroutine inmat6(elrefa, fapg, mganos)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/elrfno.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfvf.h"
#include "asterfort/mgauss.h"
!
    character(len=8), intent(in) :: elrefa, fapg
    real(kind=8), intent(out) :: mganos(MT_NBPGMX, MT_NNOMAX)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Compute the Gauss passage matrix at vertex nodes
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefa           : name of geometric support for finite element
! In  fapg             : name of Gauss integration scheme
! Out mganos           : Gauss passage matrix at vertex nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, nnos
    integer(kind=8) :: i, kp, kdim, ln, j, lm, npg, iret
    real(kind=8) :: ff(MT_NNOMAX), m(MT_NBPGMX*MT_NNOMAX)
    real(kind=8) :: p(MT_NBPGMX*MT_NNOMAX)
    real(kind=8) :: xpg(3*MT_NBPGMX), poipg(MT_NBPGMX), xg(3), det
    character(len=8) :: elref2
!
! --------------------------------------------------------------------------------------------------
!
    call elrfno(elrefa, nno, nnos, ndim)
    call elraga(elrefa, fapg, ndim, npg, xpg, poipg)
    ASSERT(npg .le. MT_NBPGMX)
!
! - Lobatto schemes => not inversible !
!
    if (fapg .eq. 'LOB5' .or. fapg .eq. 'LOB7') then
        mganos = 0.d0
        do i = 1, nnos/2
            mganos(1, i) = 1.d0
        end do
        do i = nnos/2+1, nnos
            mganos(5, i) = 1.d0
        end do
        elref2 = elrefa
        goto 100
    end if
!
! - QU4/FIS2 NON INVERSIBLE => not inversible !
!
    if (elrefa .eq. 'QU4' .and. fapg .eq. 'FIS2') then
        mganos = 0.d0
        mganos(1, 1) = 1.d0
        mganos(1, 4) = 1.d0
        mganos(2, 2) = 1.d0
        mganos(2, 3) = 1.d0
        goto 100
    end if
!
! - Get linear support
!
    if ((elrefa .eq. 'H20') .or. (elrefa .eq. 'H27')) then
        elref2 = 'HE8'
    else if ((elrefa .eq. 'P15') .or. (elrefa .eq. 'P18') .or. (elrefa .eq. 'P21')) then
        elref2 = 'PE6'
    else if ((elrefa .eq. 'HE9')) then
        elref2 = 'HE8'
    else if ((elrefa .eq. 'PE7')) then
        elref2 = 'PE6'
    else if ((elrefa .eq. 'P13') .or. (elrefa .eq. 'P19')) then
        elref2 = 'PY5'
    else if ((elrefa .eq. 'T10') .or. (elrefa .eq. 'T15')) then
        elref2 = 'TE4'
    else if ((elrefa .eq. 'TR6') .or. (elrefa .eq. 'TR7')) then
        elref2 = 'TR3'
    else if ((elrefa .eq. 'QU8') .or. (elrefa .eq. 'QU9')) then
        elref2 = 'QU4'
    else if ((elrefa .eq. 'SE3') .or. (elrefa .eq. 'SE4')) then
        elref2 = 'SE2'
    else
        elref2 = elrefa
    end if
!
!
!     CALCUL DES MATRICES M ET P :
!     ----------------------------
    do i = 1, nnos*nnos
        m(i) = 0.d0
    end do
!
    do kp = 1, npg
        do kdim = 1, ndim
            xg(kdim) = xpg(ndim*(kp-1)+kdim)
        end do
        call elrfvf(elref2, xg, ff)
        ln = (kp-1)*nnos
        do i = 1, nnos
            p(ln+i) = ff(i)
            do j = 1, nnos
                lm = nnos*(i-1)+j
                m(lm) = m(lm)+ff(i)*ff(j)
            end do
        end do
    end do
!
!     CALCUL DE LA MATRICE M-1*P :
!     ----------------------------
    call mgauss('NFVP', m, p, nnos, nnos, &
                npg, det, iret)
!
    do i = 1, nnos
        do kp = 1, npg
            mganos(kp, i) = p((kp-1)*nnos+i)
        end do
    end do
!
100 continue
!
end subroutine
