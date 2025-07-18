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
subroutine dxtloc(flex, memb, mefl, ctor, matloc)
    implicit none
#include "jeveux.h"
    real(kind=8) :: flex(*), memb(*), mefl(*), ctor
    real(kind=8) :: matloc(*)
!-----------------------------------------------------
!     IN  FLEX   : MATRICE DE FLEXION CARREE
!     IN  MEMB   : MATRICE DE MEMBRANE CARREE
!     IN  MEFL   : MATRICE MEMBRANE - FLEXION CARREE
!     IN  CTOR   : COEFF DE TORSION
!     OUT MATLOC : MATRICE DE RIGIDITE OU DE MASSE LOCALE
!                  REMPLISSAGE DE MATELEM LOCAL (171 TERMES) AVEC
!                      21 TERMES DE MEMBRANE DX DY
!                      45 TERMES DE FLEXION  DZ DRX DRY
!                      54 TERMES DE MEMBRANE/FLEXION
!                       3 TERMES DE ROTATION DRZ
!-----------------------
    integer(kind=8) :: if(45), jf(45)
    integer(kind=8) :: im(21), jm(21)
    integer(kind=8) :: ifm(36), jfm(36)
    integer(kind=8) :: imf(18), jmf(18)
    integer(kind=8) :: jz(3)
    real(kind=8) :: coef
    real(kind=8) :: cf(45), cfm(36), cmf(18)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, k
!-----------------------------------------------------------------------
    data cf/1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 1.d0, 2*-1.d0,&
     &             2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &                 1.d0, 2*-1.d0, 1.d0, 2*-1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &                 1.d0/
    data cfm/2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0,&
     &             2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0,&
     &             2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0,&
     &             2*1.d0, 2*-1.d0, 2*1.d0/
    data cmf/1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &             2*1.d0, -1.d0, 1.d0/
!     ------------------------------------------------------------------
    data jf/&
     &    6, 9, 10, 13, 14, 15, 39, 40, 41, 45,&
     &   48, 49, 50, 54, 55, 58, 59, 60, 64, 65,&
     &   66, 108, 109, 110, 114, 115, 116, 120, 123, 124,&
     &  125, 129, 130, 131, 135, 136, 139, 140, 141, 145,&
     &  146, 147, 151, 152, 153/
    data if/&
     &    1, 19, 21, 10, 12, 11, 28, 30, 29, 31,&
     &   46, 48, 47, 49, 51, 37, 39, 38, 40, 42,&
     &   41, 55, 57, 56, 58, 60, 59, 61, 73, 75,&
     &   74, 76, 78, 77, 79, 81, 64, 66, 65, 67,&
     &   69, 68, 70, 72, 71/
!     ------------------------------------------------------------------
    data jm/&
     &    1, 2, 3, 22, 23, 28, 29, 30, 35, 36,&
     &   79, 80, 85, 86, 91, 92, 93, 98, 99, 104,&
     &  105/
    data im/&
     &    1, 7, 8, 13, 14, 15, 19, 20, 21, 22,&
     &   25, 26, 27, 28, 29, 31, 32, 33, 34, 35,&
     &   36/
!     ------------------------------------------------------------------
    data jfm/&
     &    4, 5, 7, 8, 11, 12, 37, 38, 46, 47,&
     &   56, 57, 43, 44, 52, 53, 62, 63, 106, 107,&
     &  121, 122, 137, 138, 112, 113, 127, 128, 143, 144,&
     &  118, 119, 133, 134, 149, 150/
    data ifm/&
     &    1, 2, 13, 14, 7, 8, 19, 20, 31, 32,&
     &   25, 26, 21, 22, 33, 34, 27, 28, 37, 38,&
     &   49, 50, 43, 44, 39, 40, 51, 52, 45, 46,&
     &   41, 42, 53, 54, 47, 48/
!     ------------------------------------------------------------------
    data jmf/&
     &   24, 25, 26, 31, 32, 33, 81, 82, 83, 94,&
     &   95, 96, 87, 88, 89, 100, 101, 102/
    data imf/&
     &    3, 15, 9, 4, 16, 10, 5, 17, 11, 6,&
     &   18, 12, 23, 35, 29, 24, 36, 30/
!     ------------------------------------------------------------------
    data jz/21, 78, 171/
!     ------------------------------------------------------------------
!                          ---- RAZ MATLOC
    do i = 1, 171
        matloc(i) = 0.0d0
    end do
!                          ---- TERMES DE FLEXION
    do k = 1, 45
        matloc(jf(k)) = cf(k)*flex(if(k))
    end do
!                          ---- TERMES DE MEMBRANE
    do k = 1, 21
        matloc(jm(k)) = memb(im(k))
    end do
!                          ---- TERMES DE COUPLAGE FLEXION/MEMBRANE
    do k = 1, 36
        matloc(jfm(k)) = cfm(k)*mefl(ifm(k))
    end do
!                          ---- TERMES DE COUPLAGE MEMBRANE/FLEXION
    do k = 1, 18
        matloc(jmf(k)) = cmf(k)*mefl(imf(k))
    end do
!                          ---- TERMES DE ROTATION / Z
    coef = ctor*min(flex(11), flex(21), flex(41), flex(51), flex(71), flex(81))
    matloc(jz(1)) = coef
    matloc(jz(2)) = coef
    matloc(jz(3)) = coef
end subroutine
