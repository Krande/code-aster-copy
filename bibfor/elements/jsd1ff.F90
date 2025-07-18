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
subroutine jsd1ff(ip, xl, phiy, phiz, b)
    implicit none
    integer(kind=8) :: ip
    real(kind=8) :: xl, a, k, b(7, 14)
    real(kind=8) :: phiy, phiz
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES FONCTIONS DE FORME DE DEFORMATIONS
!               GENERALISEES POUTRE 7 DDL A TROIS  POINTS DE GAUSS
!
!    - ARGUMENTS:
!        DONNEES:           IP      -->   POINTS DE GAUSS
!                           XL      -->   LONGUEUR DE L'ELEMENT
!                          PHIY     -->  COEFF DE CISAILLEMENT SUIVANT Y
!                          PHIZ     -->  COEFF DE CISAILLEMENT SUIVANT Z
!
!        RESULTATS:
!                         B     <--  MATRICE D'INTERPOLATION
!
!     DESCRIPTION DE LA NUMEROTATION DU SEG2
!
!       +-----------------+
!       1                 2
!
!    L'ORDRE DE DDL EST  : 1   UX   -SUIVANT L'AXE DE LA POUTRE
!                          2   UY   -I
!                          3   UZ    I DANS LA SECTION
!
!                          4   TX   -.ROTATIONS SUIVANT OX,OY,OZ
!                          5   TY    .
!                          6   TZ    .
!                          7   TX'   .PARAMETRE DE GAUCHISSEMENT
!    DEFORMATIONS        : 1   UX'   .LOGITUDINALE......................
!                          2 UY'-TZ. .CISAILLEMENT OY
!                          3 UZ'+TY. .CISAILLEMENT  OZ
!                          4   TX'   .COURBURE TORSION
!                          5   TY'   .COURBURE FLEXION OY
!                          6   TZ'   .COURBURE FLEXION OZ
!                          7   TX''  .COURBURE GAUCHISSEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: dy, dz
!-----------------------------------------------------------------------
    a = 0.5d0*xl
!
    if (ip .eq. 1) then
        k = -sqrt(0.6d0)
    else if (ip .eq. 2) then
        k = 0.d0
    else if (ip .eq. 3) then
        k = sqrt(0.6d0)
    end if
    do i = 1, 7
        do j = 1, 14
            b(i, j) = 0
        end do
    end do
!
!                                     TIMOCH
!   UX'
!
    b(1, 1) = -0.5d0/a
    b(1, 8) = 0.5d0/a
!
!   UY'- TZ
!
    dy = 1.d0/(1.d0+phiy)
    b(2, 2) = (3.d0*dy*k**2-2.d0-dy)/(4.d0*a)-3.d0*dy*(k**2-1.d0)/(4.d0*a)
    b(2, 6) = (3.d0*dy*k**2-2.d0*k-dy)/(4.d0)-(3.d0*dy*(k**2-1.d0)+2.d0*(1.d0-k))/4.d0
    b(2, 9) = -b(2, 2)
    b(2, 13) = (3.d0*dy*k**2+2.d0*k-dy)/(4.d0)-(3.d0*dy*(k**2-1.d0)+2.d0*(1.d0+k))/4.d0
!
!   UZ' + TY
!
    dz = 1.d0/(1.d0+phiz)
    b(3, 3) = (3.d0*dz*k**2-2.d0-dz)/(4.d0*a)-3.d0*dz*(k**2-1.d0)/(4.d0*a)
    b(3, 5) = -(3.d0*dz*k**2-2.d0*k-dz)/(4.d0)+(3.d0*dz*(k**2-1.d0)+2.d0*(1.d0-k))/4.d0
    b(3, 10) = -b(3, 3)
    b(3, 12) = -(3.d0*dz*k**2+2.d0*k-dz)/(4.d0)+(3.d0*dz*(k**2-1.d0)+2.d0*(1.d0+k))/4.d0
!
! TX'
!
    b(4, 4) = (3.d0*k**2-3)/(4.d0*a)
    b(4, 7) = (3.d0*k**2-2.d0*k-1.d0)/(4.d0)
    b(4, 11) = (-3.d0*k**2+3.d0)/(4.d0*a)
    b(4, 14) = (3.d0*k**2+2.d0*k-1.d0)/(4.d0)
!
! TY'
!
    b(5, 3) = -(6.d0*dz*k)/(4.d0*a*a)
    b(5, 5) = (6.d0*dz*k-2.d0)/(4.d0*a)
    b(5, 10) = -b(5, 3)
    b(5, 12) = (6.d0*dz*k+2.d0)/(4.d0*a)
!
! TZ'
!
    b(6, 2) = +(6.d0*dy*k)/(4.d0*a*a)
    b(6, 6) = (6.d0*dy*k-2.d0)/(4.d0*a)
    b(6, 9) = -b(6, 2)
    b(6, 13) = (6.d0*dy*k+2.d0)/(4.d0*a)
!
! TX''
!
    b(7, 4) = (6.d0*k)/(4.d0*a*a)
    b(7, 7) = (6.d0*k-2.d0)/(4.d0*a)
    b(7, 11) = -b(7, 4)
    b(7, 14) = (6.d0*k+2.d0)/(4.d0*a)
end subroutine
