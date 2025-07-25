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
subroutine dpvpma(mod, imat, nbmat, tempd, materd, &
                  materf, matcst, ndt, ndi, nvi, &
                  indal)
    implicit none
#include "asterfort/rcvala.h"
    integer(kind=8) :: ndt, ndi, nvi, imat, nbmat
    real(kind=8) :: materd(nbmat, 2), materf(nbmat, 2), tempd
    character(len=8) :: mod
    character(len=3) :: matcst
! ====================================================================
! --- RECUPERATION DONNEES MATERIAU POUR DRUCKER PRAGER VISCOPLASTIQUE
! --- VISC_DRUC_PRAG -------------------------------------------------
! ====================================================================
! IN  : MOD    : TYPE DE MODELISATION -----------------------------
! --- : IMAT   : ADRESSE DU MATERIAU CODE -------------------------
! --- : NBMAT  : NOMBRE DE PARAMETRES MATERIAU --------------------
! --- : TEMPD  : TEMPERATURE BIDON --------------------------------
! OUT : MATERD : COEFFICIENTS MATERIAU A T ------------------------
! --- : MATERF : COEFFICIENTS MATERIAU A T+DT ---------------------
! ------------ : MATER(*,1) = CARACTERISTIQUES ELASTIQUES ---------
! ------------ : MATER(*,2) = CARACTERISTIQUES PLASTIQUES ---------
! --- : MATCST : 'OUI' --------------------------------------------
! --- : NDT    : NOMBRE TOTAL DE COMPOSANTES DU TENSEUR -----------
! --- : NDI    : NOMBRE DE COMPOSANTES DIRECTES DU TENSEUR --------
! --- : NVI    : NB DE VARIABLES INTERNES -------------------------
! --- : INDAL  : INDICATEUR SUR ALPHA
! =================================================================
    integer(kind=8) :: ii, indal, i, j
    real(kind=8) :: e, nu, mu, k
    real(kind=8) :: un, deux, trois
    integer(kind=8) :: cerr(17)
    character(len=16) :: nomc(17)
! =================================================================
    parameter(trois=3.0d0)
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
! =================================================================
! --- DEFINITION PARAMETRES MATERIAU ELAS -------------------------
! =================================================================
    nomc(1) = 'E        '
    nomc(2) = 'NU       '
    nomc(3) = 'ALPHA    '
    nomc(4) = 'PREF     '
    nomc(5) = 'A        '
    nomc(6) = 'N        '
    nomc(7) = 'P_PIC    '
    nomc(8) = 'P_ULT    '
    nomc(9) = 'ALPHA_0  '
    nomc(10) = 'ALPHA_PIC'
    nomc(11) = 'ALPHA_ULT'
    nomc(12) = 'R_0      '
    nomc(13) = 'R_PIC    '
    nomc(14) = 'R_ULT    '
    nomc(15) = 'BETA_0   '
    nomc(16) = 'BETA_PIC '
    nomc(17) = 'BETA_ULT '
!
    do i = 1, nbmat
        do j = 1, 2
            materd(i, j) = 0.d0
        end do
    end do
!
! =================================================================
! --- RECUPERATION DES PARAMETRES MATERIAU ------------------------
! =================================================================
    call rcvala(imat, ' ', 'ELAS', 1, 'TEMP', &
                [tempd], 3, nomc(1), materd(1, 1), cerr(1), &
                0)
    indal = 1
    if (cerr(3) .ne. 0) indal = 0
!
    call rcvala(imat, ' ', 'VISC_DRUC_PRAG', 1, 'TEMP', &
                [tempd], 14, nomc(4), materd(1, 2), cerr(4), &
                0)
! =================================================================
! - CALCUL DES MODULES DE CISAILLEMENT ET DE DEFORMATION VOLUMIQUE-
! =================================================================
    e = materd(1, 1)
    nu = materd(2, 1)
    mu = e/(deux*(un+nu))
    k = e/(trois*(un-deux*nu))
! =================================================================
! --- STOCKAGE DES MODULES CALCULES COMME PARAMETRES MATERIAU -----
! =================================================================
    materd(4, 1) = mu
    materd(5, 1) = k
! =================================================================
! --- DEFINITION D'UN MATERIAU FINAL ------------------------------
! =================================================================
    do ii = 1, nbmat
        materf(ii, 1) = materd(ii, 1)
        materf(ii, 2) = materd(ii, 2)
    end do
    matcst = 'OUI'
! =================================================================
! --- NOMBRE DE COMPOSANTES ---------------------------------------
! =================================================================
    if (mod(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
    else if ((mod(1:6) .eq. 'D_PLAN') .or. (mod(1:4) .eq. 'AXIS')) then
        ndt = 4
        ndi = 3
    end if
! =================================================================
! --- NOMBRE DE VARIABLES INTERNES --------------------------------
! =================================================================
    nvi = 4
! =================================================================
end subroutine
