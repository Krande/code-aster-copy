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
subroutine intrpo(r, s, t, nno, vh)
!
!
    implicit none
! ......................................................................
!     CALCUL DES VALEURS DES FONCTONS D'INTERPOLATION
!     ET DE LEURS DERIVEES AU POINT (R,S,T)
!
! IN    NNO  : NOMBRE DE NOEUDS
! IN    R,S,T: COORDONNES DU POINT CONSIDERE
! OUT    VH  : VALEUR DE LA FONCTION D INTERPOLATION AU POINT
! ......................................................................
!
    real(kind=8) :: vh(50), hr(3), hs(3), ht(3)
    integer(kind=8) :: n27(27)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, ih, j, k, nno
    real(kind=8) :: de, hu, r, rm, rp, s, sm
    real(kind=8) :: sp, t, tm, tp, uh, un
!-----------------------------------------------------------------------
    data un, de, hu/1.d0, 2.d0, 8.d0/
    data n27/1, 14, 5, 12, 21, 26, 4, 20, 8, 9, 15, 23, 13, 22, 27, 11, 19, 25, 2, 16,&
     &         6, 10, 17, 24, 3, 18, 7/
!
!----------HEXAEDRE A 8 NOEUDS
!
    if (nno .eq. 8) then
!
        rm = un-r
        rp = un+r
        sm = un-s
        sp = un+s
        tm = un-t
        tp = un+t
        uh = un/hu
!---------------FONCTIONS D INTERPOLATION
        vh(1) = uh*rm*sm*tm
        vh(2) = uh*rp*sm*tm
        vh(3) = uh*rp*sp*tm
        vh(4) = uh*rm*sp*tm
        vh(5) = uh*rm*sm*tp
        vh(6) = uh*rp*sm*tp
        vh(7) = uh*rp*sp*tp
        vh(8) = uh*rm*sp*tp
    end if
!
!----------ELEMENTS A 20 OU A 27 NOEUDS
!
    if (nno .eq. 20 .or. nno .eq. 27) then
        hr(1) = r*(r-un)/de
        hr(2) = un-r*r
        hr(3) = r*(r+un)/de
        hs(1) = s*(s-un)/de
        hs(2) = un-s*s
        hs(3) = s*(s+un)/de
        ht(1) = t*(t-un)/de
        ht(2) = un-t*t
        ht(3) = t*(t+un)/de
!
        i1 = 0
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    i1 = i1+1
                    ih = n27(i1)
                    vh(ih) = hr(i)*hs(j)*ht(k)
                end do
            end do
        end do
    end if
!
end subroutine
