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
subroutine dxtloc2(flex, memb, mefl, ctor, matloc)
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
!                  REMPLISSAGE DE MATELEM LOCAL (324 TERMES) AVEC
!                      36 TERMES DE MEMBRANE DX DY
!                      81 TERMES DE FLEXION  DZ DRX DRY
!                     108 TERMES DE MEMBRANE/FLEXION
!                       3 TERMES DE ROTATION DRZ
!-----------------------
    integer(kind=8) :: jf(81)
    integer(kind=8) :: jm(36)
    integer(kind=8) :: jfm(54)
    integer(kind=8) :: jmf(54)
    integer(kind=8) :: jz(3)
    real(kind=8) :: coef
    real(kind=8) :: cf(81), cfm(54)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, k
!-----------------------------------------------------------------------
    data cf/&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*-1.d0, 1.d0, 2*-1.d0, 1.d0, 2*-1.d0, 1.d0,&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*-1.d0, 1.d0, 2*-1.d0, 1.d0, 2*-1.d0, 1.d0,&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*1.d0, -1.d0, 2*1.d0, -1.d0, 2*1.d0, -1.d0,&
     &   2*-1.d0, 1.d0, 2*-1.d0, 1.d0, 2*-1.d0, 1.d0/
!     ------------------------------------------------------------------
    data cfm/12*1.d0, 6*-1.d0, 12*1.d0, 6*-1.d0, 12*1.d0, 6*-1.d0/
!     ------------------------------------------------------------------
    data jf/&
     &   39, 41, 40, 45, 47, 46, 51, 53, 52,&
     &   75, 77, 76, 81, 83, 82, 87, 89, 88,&
     &   57, 59, 58, 63, 65, 64, 69, 71, 70,&
     &  147, 149, 148, 153, 155, 154, 159, 161, 160,&
     &  183, 185, 184, 189, 191, 190, 195, 197, 196,&
     &  165, 167, 166, 171, 173, 172, 177, 179, 178,&
     &  255, 257, 256, 261, 263, 262, 267, 269, 268,&
     &  291, 293, 292, 297, 299, 298, 303, 305, 304,&
     &  273, 275, 274, 279, 281, 280, 285, 287, 286/

!     ------------------------------------------------------------------
    data jm/&
     &    1, 2, 7, 8, 13, 14,&
     &   19, 20, 25, 26, 31, 32,&
     &  109, 110, 115, 116, 121, 122,&
     &  127, 128, 133, 134, 139, 140,&
     &  217, 218, 223, 224, 229, 230,&
     &  235, 236, 241, 242, 247, 248/
!     ------------------------------------------------------------------
    data jfm/&
     &   37, 38, 43, 44, 49, 50,&
     &   73, 74, 79, 80, 85, 86,&
     &   55, 56, 61, 62, 67, 68,&
     &  145, 146, 151, 152, 157, 158,&
     &  181, 182, 187, 188, 193, 194,&
     &  163, 164, 169, 170, 175, 176,&
     &  253, 254, 259, 260, 265, 266,&
     &  289, 290, 295, 296, 301, 302,&
     &  271, 272, 277, 278, 283, 284/
!     ------------------------------------------------------------------
    data jmf/&
     &    3, 21, 111, 129, 219, 237,&
     &    5, 23, 113, 131, 221, 239,&
     &    4, 22, 112, 130, 220, 238,&
     &    9, 27, 117, 135, 225, 243,&
     &   11, 29, 119, 137, 227, 245,&
     &   10, 28, 118, 136, 226, 244,&
     &   15, 33, 123, 141, 231, 249,&
     &   17, 35, 125, 143, 233, 251,&
     &   16, 34, 124, 142, 232, 250/
!     ------------------------------------------------------------------
    data jz/96, 210, 324/
!     ------------------------------------------------------------------
!                          ---- RAZ MATLOC
    do i = 1, 324
        matloc(i) = 0.0d0
    end do
!                          ---- TERMES DE FLEXION
    do k = 1, 81
        matloc(jf(k)) = cf(k)*flex(k)
    end do
!                          ---- TERMES DE MEMBRANE
    do k = 1, 36
        matloc(jm(k)) = memb(k)
    end do
!                          ---- TERMES DE COUPLAGE FLEXION/MEMBRANE et MEMBRANE/FLEXION
    do k = 1, 54
        matloc(jfm(k)) = cfm(k)*mefl(k)
        matloc(jmf(k)) = cfm(k)*mefl(k)
    end do
!                          ---- TERMES DE ROTATION / Z
    coef = ctor*min(flex(11), flex(21), flex(41), flex(51), flex(71), flex(81))
    matloc(jz(1)) = coef
    matloc(jz(2)) = coef
    matloc(jz(3)) = coef
end subroutine
