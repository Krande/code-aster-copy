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
subroutine tranlg(nb1, nddlx, nddlet, plg, matloc, &
                  xr)
    implicit none
!
    integer(kind=8) :: nb1, nddlet, nddlx
    real(kind=8) :: xr(*)
    real(kind=8) :: matloc(nddlx, nddlx), plg(9, 3, 3)
    real(kind=8) :: matxp(48, 48), matx(51, 51)
    real(kind=8) :: kb12pt(48, 3), kb21pg(3, 48), kb22pt(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, ib, j, j1, j2
    integer(kind=8) :: jb, k, k1, kompt, nb2, nddle
!
!-----------------------------------------------------------------------
    nddle = 6*nb1
!
!     CONSTRUCTION DE LA MATRICE GLOBALE K = TLAMBDA * KB * LAMBDA
!
!             KB11 KB12                 P   0
!       KB =                  LAMBDA =
!             KB21 KB22                 0   P9
!
!       KB11(NDDLE,NDDLE) , KB12(NDDLE,3) , KB21(3,NDDLE) , KB22(3,3)
!          P(NDDLE,NDDLE) ,   P9(3,3)
!
!          P = PLG (IB,3,3) IB=1,NB1
!
!     CALCULS DU BLOC PRINCIPAL TP * KB11 * P
!
!     CONSTRUCTION DE LA MATRICE MATXP = MATLOC *  PLG  : (NDDLE,NDDLE)
!
    do i1 = 1, nddle
        do jb = 1, nb1
            do j = 1, 3
                j1 = 3*(2*jb-2)+j
                matxp(i1, j1) = matloc(i1, j1)
!
                j2 = 3*(2*jb-1)+j
                matxp(i1, j2) = 0.d0
                do k = 1, 3
                    k1 = 3*(2*jb-1)+k
                    matxp(i1, j2) = matxp(i1, j2)+matloc(i1, k1)*plg(jb, k, &
                                                                     j)
                end do
            end do
        end do
    end do
!
!   CONSTRUCTION DE LA MATRICE MATX = PLGT * MATLOC * PLG (NDDLET,NDDLE)
!
!   CONSTRUCTION EN PREMIER DE PT * KB11 * P (NDDLE,NDDLE)
!
    do ib = 1, nb1
        do i = 1, 3
            i1 = 3*(2*ib-2)+i
!
            i2 = 3*(2*ib-1)+i
            do j1 = 1, nddle
                matx(i1, j1) = matxp(i1, j1)
!
                matx(i2, j1) = 0.d0
                do k = 1, 3
                    k1 = 3*(2*ib-1)+k
!CC      MATX(I2,J1)=MATX(I2,J1)+PLGT(IB,I,K)*MATXP(K1,J1)
                    matx(i2, j1) = matx(i2, j1)+plg(ib, k, i)*matxp(k1, j1)
                end do
            end do
        end do
    end do
!
    nb2 = nb1+1
!
!     CALCULS DES SOUS MATRICES POUR FORMER LA MATRICE COMPLETE K
!
!     CALCULS DE KB * LAMBDA
!
    do i = 1, nddle
        do j = 1, 3
            kb12pt(i, j) = 0.d0
            do k = 1, 3
                k1 = nddle+k
                kb12pt(i, j) = kb12pt(i, j)+matloc(i, k1)*plg(nb2, k, j)
            end do
        end do
    end do
!
    do i = 1, 3
        i1 = nddle+i
        do jb = 1, nb1
            do j = 1, 3
                j1 = 3*(2*jb-2)+j
                kb21pg(i, j1) = matloc(i1, j1)
!
                j2 = 3*(2*jb-1)+j
                kb21pg(i, j2) = 0.d0
                do k = 1, 3
                    k1 = 3*(2*jb-1)+k
                    kb21pg(i, j2) = kb21pg(i, j2)+matloc(i1, k1)*plg(jb, k, &
                                                                     j)
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        i1 = nddle+i
        do j = 1, 3
            kb22pt(i, j) = 0.d0
            do k = 1, 3
                k1 = nddle+k
                kb22pt(i, j) = kb22pt(i, j)+matloc(i1, k1)*plg(nb2, k, j)
            end do
        end do
    end do
!
!     CALCULS DE K = TLAMBDA * KB * LAMBDA
!
    do ib = 1, nb1
        do i = 1, 3
            i1 = 3*(2*ib-2)+i
!
            i2 = 3*(2*ib-1)+i
            do j = 1, 3
                j1 = nddle+j
                matx(i1, j1) = kb12pt(i1, j)
!
                matx(i2, j1) = 0.d0
                do k = 1, 3
                    k1 = 3*(2*ib-1)+k
                    matx(i2, j1) = matx(i2, j1)+plg(ib, k, i)*kb12pt(k1, j)
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        i1 = nddle+i
        do j = 1, nddle
            matx(i1, j) = 0.d0
            do k = 1, 3
                matx(i1, j) = matx(i1, j)+plg(nb2, k, i)*kb21pg(k, j)
            end do
        end do
    end do
!
    do i = 1, 3
        i1 = nddle+i
        do j = 1, 3
            j1 = nddle+j
            matx(i1, j1) = 0.d0
            do k = 1, 3
                matx(i1, j1) = matx(i1, j1)+plg(nb2, k, i)*kb22pt(k, j)
            end do
        end do
    end do
!
!     STOCKAGE DE LA PARTIE TRIANGULAIRE SUPERIEURE DANS LE TABLEAU XR
!
    kompt = 0
    do j = 1, nddlet
        do i = 1, j
            kompt = kompt+1
            xr(kompt) = matx(i, j)
!        XR(KOMPT)=MATLOC(I,J)
        end do
    end do
!
end subroutine
