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
subroutine gddyng(kp, nno, en, x0sk, rmkm1, &
                  rmk, omkm1, ompkm1, omk, ompk, &
                  x0sec, rgmkm, rgmk, omgkm, ompgkm, &
                  omgk, ompgk)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE
!           CERTAINES GRANDEURS DYNAMIQUES AUX POINTS DE GAUSS.
!
!     IN  : KP        : NUMERO DU POINT DE GAUSS
!           NNO       : NOMBRE DE NOEUDS
!           EN        : FONCTIONS DE FORME
!           X0SK      : DERIVEE SECONDE ACTUALISEE DU DEPLACEMENT
!           RMKM1     : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I   DE L'INSTANT N+1
!           RMK       : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I+1 DE L'INSTANT N+1
!           OMKM1     : VITESSE      ANGULAIRE, A L'ITERATION I
!           OMPKM1    : ACCELERATION ANGULAIRE, A L'ITERATION I
!           OMK       : VITESSE      ANGULAIRE, A L'ITERATION I+1
!           OMPK      : ACCELERATION ANGULAIRE, A L'ITERATION I+1
!
!     OUT, AU POINT DE GAUSS NUMERO KP:
!           X0SEC     : DERIVEE SECONDE ACTUALISEE DU DEPLACEMENT
!           RGMKM     : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I   DE L'INSTANT N+1
!           RGMK      : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I+1 DE L'INSTANT N+1
!           OMGKM     : VITESSE      ANGULAIRE, A L'ITERATION I
!           OMPGKM    : ACCELERATION ANGULAIRE, A L'ITERATION I
!           OMGK      : VITESSE      ANGULAIRE, A L'ITERATION I+1
!           OMPGK     : ACCELERATION ANGULAIRE, A L'ITERATION I+1
! ------------------------------------------------------------------
    implicit none
    real(kind=8) :: en(3, 2), x0sk(3, 3), rmkm1(3, 3), rmk(3, 3), omkm1(3, 3)
    real(kind=8) :: ompkm1(3, 3), omk(3, 3), ompk(3, 3), x0sec(3), rgmkm(3)
    real(kind=8) :: rgmk(3), omgkm(3), ompgkm(3), omgk(3), ompgk(3)
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: kc, kp, ne, nno
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    zero = 0.d0
    do kc = 1, 3
        x0sec(kc) = zero
        rgmkm(kc) = zero
        rgmk(kc) = zero
        omgkm(kc) = zero
        ompgkm(kc) = zero
        omgk(kc) = zero
        ompgk(kc) = zero
    end do
    do kc = 1, 3
        do ne = 1, nno
            x0sec(kc) = x0sec(kc)+en(ne, kp)*x0sk(kc, ne)
            rgmkm(kc) = rgmkm(kc)+en(ne, kp)*rmkm1(kc, ne)
            rgmk(kc) = rgmk(kc)+en(ne, kp)*rmk(kc, ne)
            omgkm(kc) = omgkm(kc)+en(ne, kp)*omkm1(kc, ne)
            ompgkm(kc) = ompgkm(kc)+en(ne, kp)*ompkm1(kc, ne)
            omgk(kc) = omgk(kc)+en(ne, kp)*omk(kc, ne)
            ompgk(kc) = ompgk(kc)+en(ne, kp)*ompk(kc, ne)
        end do
    end do
end subroutine
