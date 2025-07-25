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
subroutine dxqloc(flex, memb, mefl, ctor, matloc)
    implicit none
#include "jeveux.h"
    real(kind=8) :: flex(*), memb(*), mefl(*), ctor
    real(kind=8) :: matloc(*)
!----------------------------------------------------------
!     IN  FLEX   : MATRICE DE FLEXION CARREE
!     IN  MEMB   : MATRICE DE MEMBRANE CARREE
!     IN  MEFL   : MATRICE MEMBRANE - FLEXION CARREE
!     IN  CTOR   : COEFF DE TORSION
!     OUT MATLOC : MATRICE DE RIGIDITE OU DE MASSE LOCALE
!                  REMPLISSAGE DE MATELEM LOCAL (300 TERMES) AVEC
!                      36 TERMES DE MEMBRANE DX DY
!                      78 TERMES DE FLEXION  DZ DRX DRY
!                      96 TERMES DE MEMBRANE/FLEXION
!                       4 TERMES DE ROTATION DRZ
!----------------------------------------------------------
    integer(kind=8) :: if(78), jf(78)
    integer(kind=8) :: im(36), jm(36)
    integer(kind=8) :: ifm(60), jfm(60)
    integer(kind=8) :: imf(36), jmf(36)
    integer(kind=8) :: jz(4)
    real(kind=8) :: coef
    real(kind=8) :: cf(78), cfm(60), cmf(36)
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
     &             2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 1.d0, 2*-1.d0,&
     &                 1.d0, 2*-1.d0, 1.d0, 2*-1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &             2*1.d0, -1.d0, 1.d0/
    data cfm/2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0,&
     &             2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0,&
     &             2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0,&
     &             2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0,&
     &             2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0, 2*1.d0,&
     &             2*-1.d0, 2*1.d0, 2*1.d0, 2*-1.d0, 2*1.d0/
    data cmf/1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &             2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0,&
     &                -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &             2*1.d0, -1.d0, 2*1.d0, -1.d0, 1.d0/
!     ------------------------------------------------------------------
    data jf/&
     &    6, 9, 10, 13, 14, 15, 39, 40, 41, 45,&
     &   48, 49, 50, 54, 55, 58, 59, 60, 64, 65,&
     &   66, 108, 109, 110, 114, 115, 116, 120, 123, 124,&
     &  125, 129, 130, 131, 135, 136, 139, 140, 141, 145,&
     &  146, 147, 151, 152, 153, 213, 214, 215, 219, 220,&
     &  221, 225, 226, 227, 231, 234, 235, 236, 240, 241,&
     &  242, 246, 247, 248, 252, 253, 256, 257, 258, 262,&
     &  263, 264, 268, 269, 270, 274, 275, 276/
    data if/&
     &    1, 25, 27, 13, 15, 14, 37, 39, 38, 40,&
     &   61, 63, 62, 64, 66, 49, 51, 50, 52, 54,&
     &   53, 73, 75, 74, 76, 78, 77, 79, 97, 99,&
     &   98, 100, 102, 101, 103, 105, 85, 87, 86, 88,&
     &   90, 89, 91, 93, 92, 109, 111, 110, 112, 114,&
     &  113, 115, 117, 116, 118, 133, 135, 134, 136, 138,&
     &  137, 139, 141, 140, 142, 144, 121, 123, 122, 124,&
     &  126, 125, 127, 129, 128, 130, 132, 131/
!     ------------------------------------------------------------------
    data jm/&
     &    1, 2, 3, 22, 23, 28, 29, 30, 35, 36,&
     &   79, 80, 85, 86, 91, 92, 93, 98, 99, 104,&
     &  105, 172, 173, 178, 179, 184, 185, 190, 191, 192,&
     &  197, 198, 203, 204, 209, 210/
    data im/&
     &    1, 9, 10, 17, 18, 19, 25, 26, 27, 28,&
     &   33, 34, 35, 36, 37, 41, 42, 43, 44, 45,&
     &   46, 49, 50, 51, 52, 53, 54, 55, 57, 58,&
     &   59, 60, 61, 62, 63, 64/
!     ------------------------------------------------------------------
    data jmf/&
     &   24, 25, 26, 31, 32, 33, 81, 82, 83, 94,&
     &   95, 96, 87, 88, 89, 100, 101, 102, 174, 175,&
     &  176, 193, 194, 195, 180, 181, 182, 199, 200, 201,&
     &  186, 187, 188, 205, 206, 207/
    data imf/&
     &    3, 19, 11, 4, 20, 12, 5, 21, 13, 6,&
     &   22, 14, 29, 45, 37, 30, 46, 38, 7, 23,&
     &   15, 8, 24, 16, 31, 47, 39, 32, 48, 40,&
     &   55, 71, 63, 56, 72, 64/
!     ------------------------------------------------------------------
    data jfm/&
     &    4, 5, 7, 8, 11, 12, 37, 38, 46, 47,&
     &   56, 57, 43, 44, 52, 53, 62, 63, 106, 107,&
     &  121, 122, 137, 138, 112, 113, 127, 128, 143, 144,&
     &  118, 119, 133, 134, 149, 150, 211, 212, 232, 233,&
     &  254, 255, 217, 218, 238, 239, 260, 261, 223, 224,&
     &  244, 245, 266, 267, 229, 230, 250, 251, 272, 273/
    data ifm/&
     &    1, 2, 17, 18, 9, 10, 25, 26, 41, 42,&
     &   33, 34, 27, 28, 43, 44, 35, 36, 49, 50,&
     &   65, 66, 57, 58, 51, 52, 67, 68, 59, 60,&
     &   53, 54, 69, 70, 61, 62, 73, 74, 89, 90,&
     &   81, 82, 75, 76, 91, 92, 83, 84, 77, 78,&
     &   93, 94, 85, 86, 79, 80, 95, 96, 87, 88/
!     ------------------------------------------------------------------
    data jz/21, 78, 171, 300/
!     ------------------------------------------------------------------
!                          ---- RAZ MATLOC
    do i = 1, 300
        matloc(i) = 0.0d0
    end do
!                          ---- TERMES DE FLEXION
    do k = 1, 78
        matloc(jf(k)) = cf(k)*flex(if(k))
    end do
!                          ---- TERMES DE MEMBRANE
    do k = 1, 36
        matloc(jm(k)) = memb(im(k))
    end do
!                          ---- TERMES DE COUPLAGE FLEXION/MEMBRANE
    do k = 1, 60
        matloc(jfm(k)) = cfm(k)*mefl(ifm(k))
    end do
!                          ---- TERMES DE COUPLAGE MEMBRANE/FLEXION
    do k = 1, 36
        matloc(jmf(k)) = cmf(k)*mefl(imf(k))
    end do
!                          ---- TERMES DE ROTATION / Z
    coef = ctor*min( &
           flex(14), flex(27), flex(53), flex(66), flex(92), flex(105), flex(131), flex(144))
    matloc(jz(1)) = coef
    matloc(jz(2)) = coef
    matloc(jz(3)) = coef
    matloc(jz(4)) = coef
end subroutine
