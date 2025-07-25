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

subroutine lcexpo(a, e, mseq, nsqu)
    implicit none
    integer(kind=8), intent(in) :: mseq, nsqu
    real(kind=8), intent(in) :: a(6)
    real(kind=8), intent(out):: e(6)
! --------------------------------------------------------------------------------------------------
!  CALCUL DE L'EXPONENTIEL D'UN TENSEUR (METHODE SCALING AND SQUARING)
! --------------------------------------------------------------------------------------------------
! A      IN  TENSEUR ARGUMENT DE L'EXPONENTIEL (1:6) AVEC RAC2
! E      OUT EXPONENTIEL DE A (1:6) AVEC RAC2
! MSEQ   IN  NOMBRE DE TERMES DANS LE DEVELOPPEMENT EN SERIE ENTIERE
! NSQU   IN  NOMBRE DE NIVEAUX DE SUBDIVISION
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: k
    real(kind=8) :: b(6), t(6)
    real(kind=8), parameter :: srac2 = 1/sqrt(2.d0)
    real(kind=8), parameter, dimension(6):: kr = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
! --------------------------------------------------------------------------------------------------
!
!
! EXPONENTIEL DE LA MATRICE DE PLUS BAS NIVEAU
!
    b = kr
    e = 0
    do k = 1, mseq
        t = b/(k*2.d0**nsqu)
        b(1) = t(1)*a(1)+0.5d0*t(4)*a(4)+0.5d0*t(5)*a(5)
        b(2) = 0.5d0*t(4)*a(4)+t(2)*a(2)+0.5d0*t(6)*a(6)
        b(3) = 0.5d0*t(5)*a(5)+0.5d0*t(6)*a(6)+t(3)*a(3)
        b(4) = t(1)*a(4)+t(4)*a(2)+srac2*t(5)*a(6)
        b(5) = t(1)*a(5)+srac2*t(4)*a(6)+t(5)*a(3)
        b(6) = srac2*t(4)*a(5)+t(2)*a(6)+t(6)*a(3)
        e = e+b
    end do
!
!
! REMONTEE DES NIVEAUX
!
    do k = 1, nsqu
        t = e
        e(1) = t(1)*t(1)+0.5d0*t(4)*t(4)+0.5d0*t(5)*t(5)
        e(2) = 0.5d0*t(4)*t(4)+t(2)*t(2)+0.5d0*t(6)*t(6)
        e(3) = 0.5d0*t(5)*t(5)+0.5d0*t(6)*t(6)+t(3)*t(3)
        e(4) = t(1)*t(4)+t(4)*t(2)+srac2*t(5)*t(6)
        e(5) = t(1)*t(5)+srac2*t(4)*t(6)+t(5)*t(3)
        e(6) = srac2*t(4)*t(5)+t(2)*t(6)+t(6)*t(3)
        e = 2*t+e
    end do
    e = e+kr
!
end subroutine
