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
subroutine matrkb(nb1, ndimx, nddlx, nddlet, ktdc, &
                  alpha, rig1, coef)
    implicit none
!
#include "asterc/r8prem.h"
    integer(kind=8) :: nb1, nb2
    integer(kind=8) :: ndimx, nddlx, nddlet
!     REAL*8 KTDC(NDDLE,NDDLE),RIG1(NDDLET,NDDLET)
!     REAL*8 ALPHA,COEF
!
    real(kind=8) :: ktdc(ndimx, ndimx), rig1(nddlx, nddlx)
    real(kind=8) :: alpha, coef
!
    real(kind=8) :: rigrl(2, 2)
    real(kind=8) :: xmin
!
    integer(kind=8) :: i1, i2, i, ib, ir, in, ii
    integer(kind=8) :: j1, j2, j, jb, jj, jr
    integer(kind=8) :: kompti, komptj
!
!
!
!     RECHERCHE DU MINIMUM DE LA MATRICE KTDC = INF DES TERMES DIAGONAUX
!
!
    xmin = 1.d0/r8prem()
!
    nb2 = nb1+1
!
    do in = 1, nb2
!
!------- ON CONSTRUIT RIGRL
!
        if (in .le. nb1) then
!
!----------- NOEUDS DE SERENDIP
            do jj = 1, 2
                jr = 5*(in-1)+jj+3
                do ii = 1, 2
                    ir = 5*(in-1)+ii+3
                    rigrl(ii, jj) = ktdc(ir, jr)
                end do
            end do
!
        else
!
!----------- SUPERNOEUD
            do jj = 1, 2
                jr = 5*nb1+jj
                do ii = 1, 2
                    ir = 5*nb1+ii
                    rigrl(ii, jj) = ktdc(ir, jr)
                end do
            end do
!
        end if
!
!-------    ON COMPARE LES DEUX PREMIERS TERMES DIAGONAUX DE RIGRL
!
        if (rigrl(1, 1) .lt. xmin) xmin = rigrl(1, 1)
        if (rigrl(2, 2) .lt. xmin) xmin = rigrl(2, 2)
!
    end do
!
!
    coef = alpha*xmin
!
!     RAIDEUR ASSOCIEE A LA ROTATION FICTIVE = COEF = ALPHA * INF
!
    do i = 1, nddlet
        do j = 1, nddlet
            rig1(i, j) = 0.d0
        end do
    end do
!
!     CONSTRUCTION DE KBARRE = KTILD EXTENDU  :  (NDDLET,NDDLET)
!
    kompti = -1
    komptj = -1
!
    nb2 = nb1+1
!
    do ib = 1, nb2
        kompti = kompti+1
        do jb = 1, nb2
            komptj = komptj+1
            if ((ib .le. nb1) .and. (jb .le. nb1)) then
                do i = 1, 5
                    i1 = 5*(ib-1)+i
                    i2 = i1+kompti
                    do j = 1, 5
                        j1 = 5*(jb-1)+j
                        j2 = j1+komptj
                        rig1(i2, j2) = ktdc(i1, j1)
                    end do
                end do
!
            else if ((ib .le. nb1) .and. (jb .eq. nb2)) then
                do i = 1, 5
                    i1 = 5*(ib-1)+i
                    i2 = i1+kompti
                    do j = 1, 2
                        j1 = 5*nb1+j
                        j2 = j1+komptj
                        rig1(i2, j2) = ktdc(i1, j1)
                    end do
                end do
!
            else if ((ib .eq. nb2) .and. (jb .le. nb1)) then
                do i = 1, 2
                    i1 = 5*nb1+i
                    i2 = i1+kompti
                    do j = 1, 5
                        j1 = 5*(jb-1)+j
                        j2 = j1+komptj
                        rig1(i2, j2) = ktdc(i1, j1)
                    end do
                end do
!
            else if ((ib .eq. nb2) .and. (jb .eq. nb2)) then
                do i = 1, 2
                    i1 = 5*nb1+i
                    i2 = i1+kompti
                    do j = 1, 2
                        j1 = 5*nb1+j
                        j2 = j1+komptj
                        rig1(i2, j2) = ktdc(i1, j1)
                    end do
                end do
!
            end if
!
        end do
!
        if (ib .le. nb1) then
            rig1(6*ib, 6*ib) = coef
        else
            rig1(nddlet, nddlet) = coef
        end if
!
        komptj = -1
    end do
!
end subroutine
