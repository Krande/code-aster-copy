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
subroutine calkbb(nno, ndim, w, def, dsidep, &
                  kbb)
! person_in_charge: sebastien.fayolle at edf.fr
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim, nno
    real(kind=8) :: w, def(2*ndim, nno, ndim)
    real(kind=8) :: kbb(ndim, ndim), dsidep(2*ndim, 2*ndim)
!-----------------------------------------------------------------------
!     BUT:  CALCUL DE LA MATRICE DE RAIDEUR LIEE A LA BULLE KBB
!-----------------------------------------------------------------------
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  W      : POIDS DU POINT DE GAUSS
! IN  DEF    : MATRICE B
! IN  DSIDEP : MATRICE TANGENTE COHERENTE POUR LA PARTIE BULLE
! OUT KBB    : MATRICE KBB
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ia, ja, na, kl, pq
    real(kind=8) :: t1
    real(kind=8) :: pbulle
    real(kind=8) :: devd(2*ndim, 2*ndim)
    real(kind=8) :: dddev(2*ndim, 2*ndim)
    real(kind=8) :: idev(6, 6), idev2(4, 4)
!
    data idev2/2.d0, -1.d0, -1.d0, 0.d0,&
     &                   -1.d0, 2.d0, -1.d0, 0.d0,&
     &                   -1.d0, -1.d0, 2.d0, 0.d0,&
     &                    0.d0, 0.d0, 0.d0, 3.d0/
    data idev/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                   -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                   -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0,&
     &                    0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0,&
     &                    0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0,&
     &                    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/
!-----------------------------------------------------------------------
!
! - INITIALISATION
    call r8inir(ndim*ndim, 0.d0, kbb, 1)
    devd(:, :) = 0.d0
    dddev(:, :) = 0.d0
!
    if (ndim .eq. 3) then
        pbulle = 4.d0
        devd(1:6, 1:6) = matmul(idev/3.d0, dsidep(1:6, 1:6))
        dddev(1:6, 1:6) = matmul(devd(1:6, 1:6), idev/3.d0)
    else if (ndim .eq. 2) then
        pbulle = 3.d0
        devd(1:4, 1:4) = matmul(idev2/3.d0, dsidep(1:4, 1:4))
        dddev(1:4, 1:4) = matmul(devd(1:4, 1:4), idev2/3.d0)
    else
        ASSERT(.false.)
    end if
!
! - CALCUL DE LA MATRICE KBB
! - BOUCLE SUR LES SOUS ELEMENTS
    do na = 1, nno
        do ia = 1, ndim
            do ja = 1, ndim
                t1 = 0.d0
                do kl = 1, 2*ndim
                    do pq = 1, 2*ndim
                        t1 = t1+def(kl, na, ia)*dddev(kl, pq)*def(pq, na, ja)
                    end do
                end do
                kbb(ia, ja) = kbb(ia, ja)+pbulle*w*t1
            end do
        end do
    end do
!
end subroutine
