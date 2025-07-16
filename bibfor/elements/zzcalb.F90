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
subroutine zzcalb(igr, iel, npg, nno, wi, &
                  desc, sig, x, y, xmin, &
                  xmax, ymin, ymax, f)
    implicit none
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, iadiel, ideb, iel, igr, ino
    integer(kind=8) :: ipg, j, nno, npg
    real(kind=8) :: xmax, xmin, ymax, ymin
!-----------------------------------------------------------------------
!
!    ESTIMATEUR ZZ (2-EME VERSION 92)
!
!                                  T
! CETTE ROUTINE CALCULE  F = SIGMA( P(X   ,Y   ) * SIG (X   ,Y   ) )
!                         I   IPG      IPG  IPG       I  IPG  IPG
!
!    AVEC P(X,Y) = (1,X,Y,XY)
!       EN PRENANT LES NNO PREMIERS MONOMES
!
!     X(IPG),Y(IPG) SONT LES COORDONNEES DES PTS DE GAUSS SUR
!     L'ELEMENT COURANT IMA
!
    real(kind=8) :: wi(1), x(1), y(1), f(9, 4), xx, yy, sig(1), b(9)
    integer(kind=8) :: desc(1)
    do ipg = 1, npg
        xx = 0.d0
        yy = 0.d0
        do ino = 1, nno
            xx = xx+wi(nno*(ipg-1)+ino)*x(ino)
            yy = yy+wi(nno*(ipg-1)+ino)*y(ino)
        end do
        xx = -1.d0+2.d0*(xx-xmin)/(xmax-xmin)
        yy = -1.d0+2.d0*(yy-ymin)/(ymax-ymin)
        b(1) = 1.d0
        b(2) = xx
        b(3) = yy
        b(4) = xx*yy
        b(5) = xx*xx
        b(6) = yy*yy
        b(7) = xx*b(4)
        b(8) = yy*b(4)
        b(9) = b(6)*b(5)
!
        ideb = desc(4+igr)
        iadiel = desc(ideb+8+4*(iel-1))
        iad = iadiel+4*(ipg-1)-1
        do i = 1, nno
            do j = 1, 4
                f(i, j) = f(i, j)+b(i)*sig(iad+j)
            end do
        end do
    end do
!
end subroutine
