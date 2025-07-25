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
subroutine lcmmsg(nomfam, nbsys, nusys, pgl2, mus, &
                  ng, mg, ir, q)
    implicit none
! person_in_charge: jean-michel.proix at edf.fr
!       IN  FAMSYS  :  NOM FAMILLE SYS GLIS
!           NUSYS   :  NUMERO SYS GLIS (FACULTATIF)
!           PGL2     :  MATRICE DE PASSAGE REPERE GLOBAL REPERE LOCAL
!           IR      :  =0 pas de rotation de reseau ; =1 : rotation
!           Q       :  matrice de rotation de reseau
!     OUT:
!           NBSYS    : NOMBRE DE SYS GLIS
!           MUS       : TENSEUR MUS POUR LE SYS GLIS NUMERO NUSYS
!
#include "asterfort/assert.h"
#include "asterfort/lcmmjs.h"
#include "asterfort/pmavec.h"
#include "asterfort/utpvlg.h"
#include "blas/dcopy.h"
    character(len=16) :: nomfam
    real(kind=8) :: mus(6), n(30, 3), m(30, 3), pgl2(3, 3), ng(3), nl(3), mg(3)
    real(kind=8) :: ml(3)
    real(kind=8) :: sqrt2, sqrt3, q(3, 3), ngr(3), mgr(3), tbsys(30, 6), norn
    real(kind=8) :: norm
    integer(kind=8) :: nbsys, nusys, k, i, j, ir
    blas_int :: b_incx, b_incy, b_n
!     ----------------------------------------------------------------
!
    if (nomfam(1:4) .eq. 'UTIL') then
        call lcmmjs(nomfam, nbsys, tbsys)
    end if
!
    if (nusys .eq. 0) then
        if (nomfam .eq. 'BCC24') then
            nbsys = 24
            goto 150
!
        else if (nomfam .eq. 'OCTAEDRIQUE') then
            nbsys = 12
            goto 150
!
        else if (nomfam .eq. 'CUBIQUE1') then
            nbsys = 12
            goto 150
!
        else if (nomfam .eq. 'CUBIQUE2') then
            nbsys = 12
            goto 150
!
        else if (nomfam .eq. 'ZIRCONIUM') then
            nbsys = 30
            goto 150
!
        else if (nomfam .eq. 'UNIAXIAL') then
            nbsys = 1
            goto 150
!
        else if (nomfam(1:4) .eq. 'UTIL') then
            goto 150
!
        else
            ASSERT(.false.)
        end if
    end if
!
    sqrt2 = sqrt(2.d0)
    sqrt3 = sqrt(3.d0)
    mus(:) = 0.d0
    if (nomfam .eq. 'ZIRCONIUM') then
!  prism1
        n(1, 1) = 1.0d0
        n(1, 2) = 0.0d0
        n(1, 3) = 0.0d0
        n(2, 1) = -0.5d0
        n(2, 2) = 0.866025403784d0
        n(2, 3) = 0.0d0
        n(3, 1) = -0.5d0
        n(3, 2) = -0.866025403784d0
        n(3, 3) = 0.0d0
        m(1, 1) = 0.0d0
        m(1, 2) = -1.0d0
        m(1, 3) = 0.0d0
        m(2, 1) = 0.866025403784d0
        m(2, 2) = 0.5d0
        m(2, 3) = 0.0d0
        m(3, 1) = -0.866025403784d0
        m(3, 2) = 0.5d0
        m(3, 3) = 0.0d0
!
!  basal1
        n(4, 1) = 1.0d0
        n(4, 2) = 0.0d0
        n(4, 3) = 0.0d0
        n(5, 1) = -0.5d0
        n(5, 2) = 0.866025403784d0
        n(5, 3) = 0.0d0
        n(6, 1) = -0.5d0
        n(6, 2) = -0.866025403784d0
        n(6, 3) = 0.0d0
        m(4, 1) = 0.0d0
        m(4, 2) = 0.0d0
        m(4, 3) = 1.0d0
        m(5, 1) = 0.0d0
        m(5, 2) = 0.0d0
        m(5, 3) = 1.0d0
        m(6, 1) = 0.0d0
        m(6, 2) = 0.0d0
        m(6, 3) = 1.0d0
!
!  pyr_a1
        n(7, 1) = 1.0d0
        n(7, 2) = 0.0d0
        n(7, 3) = 0.0d0
        n(8, 1) = 1.0d0
        n(8, 2) = 0.0d0
        n(8, 3) = 0.0d0
        n(9, 1) = -0.5d0
        n(9, 2) = 0.866025403784d0
        n(9, 3) = 0.0d0
        n(10, 1) = -0.5d0
        n(10, 2) = 0.866025403784d0
        n(10, 3) = 0.0d0
        n(11, 1) = -0.5d0
        n(11, 2) = -0.866025403784d0
        n(11, 3) = 0.0d0
        n(12, 1) = -0.5d0
        n(12, 2) = -0.866025403784d0
        n(12, 3) = 0.0d0
        m(7, 1) = 0.0d0
        m(7, 2) = -0.87856329157d0
        m(7, 3) = 0.477625944339d0
        m(8, 1) = 0.0d0
        m(8, 2) = 0.87856329157d0
        m(8, 3) = 0.477625944339d0
        m(9, 1) = -0.760858129332d0
        m(9, 2) = -0.439281645785d0
        m(9, 3) = 0.477625944339d0
        m(10, 1) = 0.760858129332d0
        m(10, 2) = 0.439281645785d0
        m(10, 3) = 0.477625944339d0
        m(11, 1) = -0.760858129332d0
        m(11, 2) = 0.439281645785d0
        m(11, 3) = 0.477625944339d0
        m(12, 1) = 0.760858129332d0
        m(12, 2) = -0.439281645785d0
        m(12, 3) = 0.477625944339d0
!
!  pyr_c_a1
        n(13, 1) = 0.531670580449d0
        n(13, 2) = 0.0d0
        n(13, 3) = -0.846951234656d0
        n(14, 1) = 0.265835290225d0
        n(14, 2) = 0.460440229114d0
        n(14, 3) = -0.846951234656d0
        n(15, 1) = 0.265835290225d0
        n(15, 2) = 0.460440229114d0
        n(15, 3) = -0.846951234656d0
        n(16, 1) = -0.265835290225d0
        n(16, 2) = 0.460440229114d0
        n(16, 3) = -0.846951234656d0
        n(17, 1) = -0.265835290225d0
        n(17, 2) = 0.460440229114d0
        n(17, 3) = -0.846951234656d0
        n(18, 1) = -0.531670580449d0
        n(18, 2) = 0.0d0
        n(18, 3) = -0.846951234656d0
        n(19, 1) = -0.531670580449d0
        n(19, 2) = 0.0d0
        n(19, 3) = -0.846951234656d0
        n(20, 1) = -0.265835290225d0
        n(20, 2) = -0.460440229114d0
        n(20, 3) = -0.846951234656d0
        n(21, 1) = -0.265835290225d0
        n(21, 2) = -0.460440229114d0
        n(21, 3) = -0.846951234656d0
        n(22, 1) = 0.265835290225d0
        n(22, 2) = -0.460440229114d0
        n(22, 3) = -0.846951234656d0
        n(23, 1) = 0.265835290225d0
        n(23, 2) = -0.460440229114d0
        n(23, 3) = -0.846951234656d0
        n(24, 1) = 0.531670580449d0
        n(24, 2) = 0.0d0
        n(24, 3) = -0.846951234656d0
        m(13, 1) = 0.760858129332d0
        m(13, 2) = 0.439281645785d0
        m(13, 3) = 0.477625944339d0
        m(14, 1) = 0.760858129332d0
        m(14, 2) = 0.439281645785d0
        m(14, 3) = 0.477625944339d0
        m(15, 1) = 0.0d0
        m(15, 2) = 0.87856329157d0
        m(15, 3) = 0.477625944339d0
        m(16, 1) = 0.0d0
        m(16, 2) = 0.87856329157d0
        m(16, 3) = 0.477625944339d0
        m(17, 1) = -0.760858129332d0
        m(17, 2) = 0.439281645785d0
        m(17, 3) = 0.477625944339d0
        m(18, 1) = -0.760858129332d0
        m(18, 2) = 0.439281645785d0
        m(18, 3) = 0.477625944339d0
        m(19, 1) = -0.760858129332d0
        m(19, 2) = -0.439281645785d0
        m(19, 3) = 0.477625944339d0
        m(20, 1) = -0.760858129332d0
        m(20, 2) = -0.439281645785d0
        m(20, 3) = 0.477625944339d0
        m(21, 1) = 0.0d0
        m(21, 2) = -0.87856329157d0
        m(21, 3) = 0.477625944339d0
        m(22, 1) = 0.0d0
        m(22, 2) = -0.87856329157d0
        m(22, 3) = 0.477625944339d0
        m(23, 1) = 0.760858129332d0
        m(23, 2) = -0.439281645785d0
        m(23, 3) = 0.477625944339d0
        m(24, 1) = 0.760858129332d0
        m(24, 2) = -0.439281645785d0
        m(24, 3) = 0.477625944339d0
!
!  pyr2_c_a1
        n(25, 1) = -0.265835290225d0
        n(25, 2) = -0.460440229114d0
        n(25, 3) = 0.846951234656d0
        n(26, 1) = 0.531670580449d0
        n(26, 2) = 0.0d0
        n(26, 3) = 0.846951234656d0
        n(27, 1) = -0.265835290225d0
        n(27, 2) = 0.460440229114d0
        n(27, 3) = 0.846951234656d0
        n(28, 1) = 0.265835290225d0
        n(28, 2) = 0.460440229114d0
        n(28, 3) = 0.846951234656d0
        n(29, 1) = -0.531670580449d0
        n(29, 2) = 0.0d0
        n(29, 3) = 0.846951234656d0
        n(30, 1) = 0.265835290225d0
        n(30, 2) = -0.460440229114d0
        n(30, 3) = 0.846951234656d0
        m(25, 1) = 0.423475617328d0
        m(25, 2) = 0.733481284978d0
        m(25, 3) = 0.531670580449d0
        m(26, 1) = -0.846951234656d0
        m(26, 2) = 0.0d0
        m(26, 3) = 0.531670580449d0
        m(27, 1) = 0.423475617328d0
        m(27, 2) = -0.733481284978d0
        m(27, 3) = 0.531670580449d0
        m(28, 1) = -0.423475617328d0
        m(28, 2) = -0.733481284978d0
        m(28, 3) = 0.531670580449d0
        m(29, 1) = 0.846951234656d0
        m(29, 2) = 0.0d0
        m(29, 3) = 0.531670580449d0
        m(30, 1) = -0.423475617328d0
        m(30, 2) = 0.733481284978d0
        m(30, 3) = 0.531670580449d0
!
!        N ET L SONT UNITAIRES
    else if (nomfam .eq. 'OCTAEDRIQUE') then
! FCC LATTICE
        n(1, 1) = 1.d0
        n(1, 2) = 1.d0
        n(1, 3) = 1.d0
        n(2, 1) = 1.d0
        n(2, 2) = 1.d0
        n(2, 3) = 1.d0
        n(3, 1) = 1.d0
        n(3, 2) = 1.d0
        n(3, 3) = 1.d0
        n(4, 1) = 1.d0
        n(4, 2) = -1.d0
        n(4, 3) = 1.d0
        n(5, 1) = 1.d0
        n(5, 2) = -1.d0
        n(5, 3) = 1.d0
        n(6, 1) = 1.d0
        n(6, 2) = -1.d0
        n(6, 3) = 1.d0
        n(7, 1) = -1.d0
        n(7, 2) = 1.d0
        n(7, 3) = 1.d0
        n(8, 1) = -1.d0
        n(8, 2) = 1.d0
        n(8, 3) = 1.d0
        n(9, 1) = -1.d0
        n(9, 2) = 1.d0
        n(9, 3) = 1.d0
        n(10, 1) = -1.d0
        n(10, 2) = -1.d0
        n(10, 3) = 1.d0
        n(11, 1) = -1.d0
        n(11, 2) = -1.d0
        n(11, 3) = 1.d0
        n(12, 1) = -1.d0
        n(12, 2) = -1.d0
        n(12, 3) = 1.d0
        m(1, 1) = -1.d0
        m(1, 2) = 0.d0
        m(1, 3) = 1.d0
        m(2, 1) = 0.d0
        m(2, 2) = -1.d0
        m(2, 3) = 1.d0
        m(3, 1) = -1.d0
        m(3, 2) = 1.d0
        m(3, 3) = 0.d0
        m(4, 1) = -1.d0
        m(4, 2) = 0.d0
        m(4, 3) = 1.d0
        m(5, 1) = 0.d0
        m(5, 2) = 1.d0
        m(5, 3) = 1.d0
        m(6, 1) = 1.d0
        m(6, 2) = 1.d0
        m(6, 3) = 0.d0
        m(7, 1) = 0.d0
        m(7, 2) = -1.d0
        m(7, 3) = 1.d0
        m(8, 1) = 1.d0
        m(8, 2) = 1.d0
        m(8, 3) = 0.d0
        m(9, 1) = 1.d0
        m(9, 2) = 0.d0
        m(9, 3) = 1.d0
        m(10, 1) = -1.d0
        m(10, 2) = 1.d0
        m(10, 3) = 0.d0
        m(11, 1) = 1.d0
        m(11, 2) = 0.d0
        m(11, 3) = 1.d0
        m(12, 1) = 0.d0
        m(12, 2) = 1.d0
        m(12, 3) = 1.d0
!        N ET L DOIVENT ETRE UNITAIRES
        do j = 1, 12
            do k = 1, 3
                m(j, k) = m(j, k)/sqrt2
                n(j, k) = n(j, k)/sqrt3
            end do
        end do
    else if (nomfam .eq. 'CUBIQUE1') then
!        BCC LATTICE, {110} SLIP
        n(1, 1) = 1.d0
        n(1, 2) = 1.d0
        n(1, 3) = 0.d0
        n(2, 1) = -1.d0
        n(2, 2) = 0.d0
        n(2, 3) = 1.d0
        n(3, 1) = 0.d0
        n(3, 2) = -1.d0
        n(3, 3) = -1.d0
        n(4, 1) = 0.d0
        n(4, 2) = -1.d0
        n(4, 3) = 1.d0
        n(5, 1) = 1.d0
        n(5, 2) = 0.d0
        n(5, 3) = -1.d0
        n(6, 1) = -1.d0
        n(6, 2) = 1.d0
        n(6, 3) = 0.d0
        n(7, 1) = -1.d0
        n(7, 2) = -1.d0
        n(7, 3) = 0.d0
        n(8, 1) = 1.d0
        n(8, 2) = 0.d0
        n(8, 3) = 1.d0
        n(9, 1) = 0.d0
        n(9, 2) = 1.d0
        n(9, 3) = -1.d0
        n(10, 1) = 1.d0
        n(10, 2) = -1.d0
        n(10, 3) = 0.d0
        n(11, 1) = -1.d0
        n(11, 2) = 0.d0
        n(11, 3) = -1.d0
        n(12, 1) = 0.d0
        n(12, 2) = 1.d0
        n(12, 3) = 1.d0
!
        m(1, 1) = 1.d0
        m(1, 2) = -1.d0
        m(1, 3) = 1.d0
        m(2, 1) = 1.d0
        m(2, 2) = -1.d0
        m(2, 3) = 1.d0
        m(3, 1) = 1.d0
        m(3, 2) = -1.d0
        m(3, 3) = 1.d0
        m(4, 1) = 1.d0
        m(4, 2) = 1.d0
        m(4, 3) = 1.d0
        m(5, 1) = 1.d0
        m(5, 2) = 1.d0
        m(5, 3) = 1.d0
        m(6, 1) = 1.d0
        m(6, 2) = 1.d0
        m(6, 3) = 1.d0
        m(7, 1) = -1.d0
        m(7, 2) = 1.d0
        m(7, 3) = 1.d0
        m(8, 1) = -1.d0
        m(8, 2) = 1.d0
        m(8, 3) = 1.d0
        m(9, 1) = -1.d0
        m(9, 2) = 1.d0
        m(9, 3) = 1.d0
        m(10, 1) = 1.d0
        m(10, 2) = 1.d0
        m(10, 3) = -1.d0
        m(11, 1) = 1.d0
        m(11, 2) = 1.d0
        m(11, 3) = -1.d0
        m(12, 1) = 1.d0
        m(12, 2) = 1.d0
        m(12, 3) = -1.d0
!
!        N ET L DOIVENT ETRE UNITAIRES
        do j = 1, 12
            do k = 1, 3
                m(j, k) = m(j, k)/sqrt3
                n(j, k) = n(j, k)/sqrt2
            end do
        end do
    else if (nomfam .eq. 'CUBIQUE2') then
!        BCC LATTICE, {211} SLIP
        n(1, 1) = 2.d0
        n(1, 2) = -1.d0
        n(1, 3) = 1.d0
        n(2, 1) = 1.d0
        n(2, 2) = -2.d0
        n(2, 3) = -1.d0
        n(3, 1) = 1.d0
        n(3, 2) = 1.d0
        n(3, 3) = 2.d0
        n(4, 1) = 2.d0
        n(4, 2) = 1.d0
        n(4, 3) = 1.d0
        n(5, 1) = 1.d0
        n(5, 2) = 2.d0
        n(5, 3) = -1.d0
        n(6, 1) = 1.d0
        n(6, 2) = -1.d0
        n(6, 3) = 2.d0
        n(7, 1) = 2.d0
        n(7, 2) = 1.d0
        n(7, 3) = -1.d0
        n(8, 1) = 1.d0
        n(8, 2) = 2.d0
        n(8, 3) = 1.d0
        n(9, 1) = 1.d0
        n(9, 2) = -1.d0
        n(9, 3) = -2.d0
        n(10, 1) = 2.d0
        n(10, 2) = -1.d0
        n(10, 3) = -1.d0
        n(11, 1) = 1.d0
        n(11, 2) = -2.d0
        n(11, 3) = 1.d0
        n(12, 1) = 1.d0
        n(12, 2) = 1.d0
        n(12, 3) = -2.d0
        m(1, 1) = 1.d0
        m(1, 2) = 1.d0
        m(1, 3) = -1.d0
        m(2, 1) = 1.d0
        m(2, 2) = 1.d0
        m(2, 3) = -1.d0
        m(3, 1) = 1.d0
        m(3, 2) = 1.d0
        m(3, 3) = -1.d0
        m(4, 1) = 1.d0
        m(4, 2) = -1.d0
        m(4, 3) = -1.d0
        m(5, 1) = 1.d0
        m(5, 2) = -1.d0
        m(5, 3) = -1.d0
        m(6, 1) = 1.d0
        m(6, 2) = -1.d0
        m(6, 3) = -1.d0
        m(7, 1) = 1.d0
        m(7, 2) = -1.d0
        m(7, 3) = 1.d0
        m(8, 1) = 1.d0
        m(8, 2) = -1.d0
        m(8, 3) = 1.d0
        m(9, 1) = 1.d0
        m(9, 2) = -1.d0
        m(9, 3) = 1.d0
        m(10, 1) = 1.d0
        m(10, 2) = 1.d0
        m(10, 3) = 1.d0
        m(11, 1) = 1.d0
        m(11, 2) = 1.d0
        m(11, 3) = 1.d0
        m(12, 1) = 1.d0
        m(12, 2) = 1.d0
        m(12, 3) = 1.d0
!        N ET L DOIVENT ETRE UNITAIRES
        do j = 1, 12
            do k = 1, 3
                m(j, k) = m(j, k)/sqrt3
                n(j, k) = n(j, k)/sqrt(6.d0)
            end do
        end do
    else if (nomfam .eq. 'BCC24') then
!        BCC LATTICE, {110} SLIP
        n(1, 1) = 0.d0
        n(1, 2) = 1.d0
        n(1, 3) = 1.d0
        n(2, 1) = 1.d0
        n(2, 2) = 0.d0
        n(2, 3) = 1.d0
        n(3, 1) = 1.d0
        n(3, 2) = -1.d0
        n(3, 3) = 0.d0
        n(4, 1) = 0.d0
        n(4, 2) = 1.d0
        n(4, 3) = -1.d0
        n(5, 1) = 1.d0
        n(5, 2) = 0.d0
        n(5, 3) = 1.d0
        n(6, 1) = 1.d0
        n(6, 2) = 1.d0
        n(6, 3) = 0.d0
        n(7, 1) = 0.d0
        n(7, 2) = 1.d0
        n(7, 3) = 1.d0
        n(8, 1) = 1.d0
        n(8, 2) = 0.d0
        n(8, 3) = -1.d0
        n(9, 1) = 1.d0
        n(9, 2) = 1.d0
        n(9, 3) = 0.d0
        n(10, 1) = 0.d0
        n(10, 2) = 1.d0
        n(10, 3) = -1.d0
        n(11, 1) = 1.d0
        n(11, 2) = 0.d0
        n(11, 3) = -1.d0
        n(12, 1) = 1.d0
        n(12, 2) = -1.d0
        n(12, 3) = 0.d0
        m(1, 1) = 1.d0
        m(1, 2) = 1.d0
        m(1, 3) = -1.d0
        m(2, 1) = 1.d0
        m(2, 2) = 1.d0
        m(2, 3) = -1.d0
        m(3, 1) = 1.d0
        m(3, 2) = 1.d0
        m(3, 3) = -1.d0
        m(4, 1) = 1.d0
        m(4, 2) = -1.d0
        m(4, 3) = -1.d0
        m(5, 1) = 1.d0
        m(5, 2) = -1.d0
        m(5, 3) = -1.d0
        m(6, 1) = 1.d0
        m(6, 2) = -1.d0
        m(6, 3) = -1.d0
        m(7, 1) = 1.d0
        m(7, 2) = -1.d0
        m(7, 3) = 1.d0
        m(8, 1) = 1.d0
        m(8, 2) = -1.d0
        m(8, 3) = 1.d0
        m(9, 1) = 1.d0
        m(9, 2) = -1.d0
        m(9, 3) = 1.d0
        m(10, 1) = 1.d0
        m(10, 2) = 1.d0
        m(10, 3) = 1.d0
        m(11, 1) = 1.d0
        m(11, 2) = 1.d0
        m(11, 3) = 1.d0
        m(12, 1) = 1.d0
        m(12, 2) = 1.d0
        m(12, 3) = 1.d0
        do j = 1, 12
            do k = 1, 3
                m(j, k) = m(j, k)/sqrt3
                n(j, k) = n(j, k)/sqrt2
            end do
        end do
!        BCC LATTICE, {211} SLIP
        n(13, 1) = 2.d0
        n(13, 2) = -1.d0
        n(13, 3) = 1.d0
        n(14, 1) = 1.d0
        n(14, 2) = -2.d0
        n(14, 3) = -1.d0
        n(15, 1) = 1.d0
        n(15, 2) = 1.d0
        n(15, 3) = 2.d0
        n(16, 1) = 2.d0
        n(16, 2) = 1.d0
        n(16, 3) = 1.d0
        n(17, 1) = 1.d0
        n(17, 2) = 2.d0
        n(17, 3) = -1.d0
        n(18, 1) = 1.d0
        n(18, 2) = -1.d0
        n(18, 3) = 2.d0
        n(19, 1) = 2.d0
        n(19, 2) = 1.d0
        n(19, 3) = -1.d0
        n(20, 1) = 1.d0
        n(20, 2) = 2.d0
        n(20, 3) = 1.d0
        n(21, 1) = 1.d0
        n(21, 2) = -1.d0
        n(21, 3) = -2.d0
        n(22, 1) = 2.d0
        n(22, 2) = -1.d0
        n(22, 3) = -1.d0
        n(23, 1) = 1.d0
        n(23, 2) = -2.d0
        n(23, 3) = 1.d0
        n(24, 1) = 1.d0
        n(24, 2) = 1.d0
        n(24, 3) = -2.d0
!
        m(13, 1) = 1.d0
        m(13, 2) = 1.d0
        m(13, 3) = -1.d0
        m(14, 1) = 1.d0
        m(14, 2) = 1.d0
        m(14, 3) = -1.d0
        m(15, 1) = 1.d0
        m(15, 2) = 1.d0
        m(15, 3) = -1.d0
        m(16, 1) = 1.d0
        m(16, 2) = -1.d0
        m(16, 3) = -1.d0
        m(17, 1) = 1.d0
        m(17, 2) = -1.d0
        m(17, 3) = -1.d0
        m(18, 1) = 1.d0
        m(18, 2) = -1.d0
        m(18, 3) = -1.d0
        m(19, 1) = 1.d0
        m(19, 2) = -1.d0
        m(19, 3) = 1.d0
        m(20, 1) = 1.d0
        m(20, 2) = -1.d0
        m(20, 3) = 1.d0
        m(21, 1) = 1.d0
        m(21, 2) = -1.d0
        m(21, 3) = 1.d0
        m(22, 1) = 1.d0
        m(22, 2) = 1.d0
        m(22, 3) = 1.d0
        m(23, 1) = 1.d0
        m(23, 2) = 1.d0
        m(23, 3) = 1.d0
        m(24, 1) = 1.d0
        m(24, 2) = 1.d0
        m(24, 3) = 1.d0
!        N ET L DOIVENT ETRE UNITAIRES
        do j = 13, 24
            do k = 1, 3
                m(j, k) = m(j, k)/sqrt3
                n(j, k) = n(j, k)/sqrt(6.d0)
            end do
        end do
    else if (nomfam .eq. 'UNIAXIAL') then
        n(1, 1) = 1.d0
        n(1, 2) = 0.d0
        n(1, 3) = 0.d0
        m(1, 1) = 1.d0
        m(1, 2) = 0.d0
        m(1, 3) = 0.d0
    else if (nomfam(1:4) .eq. 'UTIL') then
        do i = 1, nbsys
            norn = sqrt(tbsys(i, 1)**2+tbsys(i, 2)**2+tbsys(i, 3)**2)
            norm = sqrt(tbsys(i, 4)**2+tbsys(i, 5)**2+tbsys(i, 6)**2)
            do j = 1, 3
                n(i, j) = tbsys(i, j)/norn
                m(i, j) = tbsys(i, j+3)/norm
            end do
        end do
!
    end if
!     POUR LE SYSTEME K, EXPRESSION DE N ET L DANS REPERE GLOBAL
    k = nusys
    do j = 1, 3
        nl(j) = n(k, j)
        ml(j) = m(k, j)
    end do
    call utpvlg(1, 3, pgl2, nl, ng)
    call utpvlg(1, 3, pgl2, ml, mg)
!     rotation de reseau
    if (ir .eq. 1) then
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ng, b_incx, ngr, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, mg, b_incx, mgr, b_incy)
        call pmavec('ZERO', 3, q, ngr, ng)
        call pmavec('ZERO', 3, q, mgr, mg)
    else
        ASSERT(ir .eq. 0)
    end if
    do j = 1, 3
        mus(j) = ng(j)*mg(j)
    end do
!     SQRT(2) PAR HOMOGENEITE AVEC NMPL3D.
    mus(4) = 0.5d0*(ng(1)*mg(2)+ng(2)*mg(1))*sqrt2
    mus(5) = 0.5d0*(ng(1)*mg(3)+ng(3)*mg(1))*sqrt2
    mus(6) = 0.5d0*(ng(2)*mg(3)+ng(3)*mg(2))*sqrt2
!
150 continue
end subroutine
