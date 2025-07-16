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
subroutine facint(nbpas, dim, longh, vec1, vec2, &
                  long, s, r, d, u, &
                  v, w)
    implicit none
!    PREPARATION DES MATRICES INTERSPECTRALES POUR FACTORISATION
! ----------------------------------------------------------------------
!          NBPAS : NOMBRE DA PAS DE DISCRETISATION DE LA MATRICE INTERSP
!          DIM   : DIMENSION DE LA MATRICE
!          DIMH  : NOMBRE DE FONCTIONS DECRIVANT LA MATRICE
!    IN  : VEC1  : VECTEUR DES VALEURS DES FONCTIONS AVANT FACTORISATION
!    OUT : VEC2  : VECTEUR DES VALEURS DES FONCTIONS APRES FACTORISATION
!          LONG  : LONGUEUR DES VECTEURS VEC1 ET VEC2
!             S  : MATRICE DE TRAVAIL DE DIMENSION DIM, A FACTORISER
!             R  : MATRICE DE TRAVAIL, RESULTAT DE LA FACTORISATION
!             D  : VECTEUR DE TRAVAIL
#include "jeveux.h"
#include "asterfort/diaghr.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: dim, long, longh
    complex(kind=8) :: s(dim, dim), r(dim, dim), u(*), w(*)
    real(kind=8) :: d(dim), vec1(long), vec2(longh), v(*)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ico, icomp, ix, iy, j, k
    integer(kind=8) :: l, nbpas, nbpt1, nbpt2
    real(kind=8) :: ai, ar, az, bz, si, sr, uu
!
!-----------------------------------------------------------------------
    nbpt1 = nbpas
    nbpt2 = nbpt1*2
    do l = 1, nbpt1
        icomp = 0
        do j = 1, dim
            do i = 1, j
                icomp = icomp+1
                ix = l+(icomp-1)*nbpt2+nbpt1
                iy = ix+nbpt1
                s(i, j) = dcmplx(vec1(ix), vec1(iy))
                if (i .ne. j) then
                    s(j, i) = dconjg(s(i, j))
                end if
            end do
        end do
        sr = dble(s(1, 1))
        si = dimag(s(1, 1))
        if (sr .eq. 0.d0 .and. si .eq. 0.d0) then
            r(1, 1) = s(1, 1)
            do i = 1, dim
                do j = 1, dim
                    sr = dble(s(i, j))
                    si = dimag(s(i, j))
                    if (sr .ne. 0.d0 .or. si .ne. 0.d0) then
                        call utmess('F', 'ALGORITH3_60')
                    end if
                    r(i, j) = s(i, j)
                end do
            end do
        else
!
!     --- FACTORISATION ---
!
            call diaghr(dim, s, dim, d, r, &
                        dim, u, v, w)
!
            do j = 1, dim
                uu = 0.d0
                do i = 1, dim
                    ar = dble(r(i, j))
                    ai = dimag(r(i, j))
                    uu = ar*ar+ai*ai+uu
                end do
                uu = sqrt(uu)
                do k = 1, dim
                    az = dble(r(k, j))/uu
                    bz = dimag(r(k, j))/uu
                    r(k, j) = dcmplx(az, bz)
                end do
            end do
            do i = 1, dim
                do j = 1, dim
                    if (d(j) .lt. 0.d0) then
                        d(j) = 0.d0
                    end if
                    r(i, j) = r(i, j)*sqrt(d(j))
                end do
            end do
        end if
        ico = 0
        do j = 1, dim
            do i = 1, dim
                ico = ico+1
                ix = l+(ico-1)*nbpt2+nbpt1
                iy = ix+nbpt1
                vec2(ix) = dble(r(i, j))
                vec2(iy) = dimag(r(i, j))
            end do
        end do
    end do
end subroutine
