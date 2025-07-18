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

subroutine mcpstr(stress, tridim, pstrs, eigprj, ii, &
                  jj, mm, codret)
! ----------------------------------------------------------------------
!
! OBJECT: COMPUTE EIGENVALUES OF STRESS AND RE-ORDER
!
! IN  STRESS  : CONTRAINTE
! IN  TRIDIM  : [LOGICAL] LA MODELISATION EST-ELLE 3D?
!
! OUT PSTRS   : CONTRAINTES PRINCIPALES (DIM=3)
! OUT EIGPRJ  : DIRECTIONS PRINCIPALES  (DIM=3xNDIM)
! OUT II      : INDICE DE LA CONTRAINTE PRINCIPALE MINEURE
! OUT JJ      : INDICE DE LA CONTRAINTE PRINCIPALE MAJEURE
! OUT MM      : INDICE DE LA CONTRAINTE PRINCIPALE INTERMEDIAIRE
! OUT CODRET  : CODE RETOUR
!               = | 0: OK
!                 | 1: NOOK
!
! ----------------------------------------------------------------------
    implicit none
! ======================================================================
!
    integer(kind=8)      :: codret
!
#include "asterf_types.h"
#include "asterfort/jacobi.h"
!
! Declaration of constant parameters
    integer(kind=8)      :: mmax, nmax, ndt, ndi
    parameter(mmax=3, nmax=6)
!
! Declaration of real type variables
    real(kind=8) :: stress(nmax), vaux(mmax)
    real(kind=8) :: pstrs1, pstrs2, pstrs3
    real(kind=8) :: eigprj(mmax, mmax), pstrs(mmax)
    real(kind=8) :: r0, r1, r2, r3, r4, small, tol
    real(kind=8) :: tu(nmax), tr(nmax), t1(nmax)
!
! Declaration of integer type variables
    integer(kind=8) :: itri, iorder, ii, jj, mm
    integer(kind=8) :: mxiter, i, itjac1
!
! Declaration of integer type variables
    aster_logical :: tridim
!
! Declaration of constant variables
    data r0, r1, r2, r3, r4, small, tol/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 1.d-06, 1.d-10/
    data mxiter/50/
!
! Declaration of Common space variables
    common/tdim/ndt, ndi
!
    codret = 0
!
! Initialize unit matrix = (1 0 0 1 0 1) for Jacobi
    t1(:) = r0
    t1(1) = r1
    t1(4) = r1
    t1(6) = r1
!
! Spectral decomposition of the trial stress
!
! ITRI =  0 : TRI EN VALEUR RELATIVE
!         1 : TRI EN VALEUR ABSOLUE
!         2 : PAS DE TRI
    itri = 2
! IORDER =  0 : TRI PAR ORDRE CROISSANT
!           1 : TRI PAR ORDRE DECROISSANT
!           2 : PAS DE TRI
    iorder = 2
! Matrix  TR = (SIXX SIXY SIXZ SIYY SIYZ SIZZ) for Jacobi
! Produce EIGPRJ: Base Projection Matrix from initial base
!                 to principal directions base
!         PSTRS : principal stresses
    if (tridim) then
        tr(1) = stress(1)
        tr(2) = stress(4)
        tr(3) = stress(5)
        tr(4) = stress(2)
        tr(5) = stress(6)
        tr(6) = stress(3)
    else
        tr(1) = stress(1)
        tr(2) = stress(4)
        tr(3) = r0
        tr(4) = stress(2)
        tr(5) = r0
        tr(6) = r0
    end if
!
! Unit matrix = (1 0 0 1 0 1) for Jacobi
    tu(1:nmax) = t1(1:nmax)
!
    eigprj(:, :) = r0
!
    call jacobi(mmax, mxiter, tol, tol, tr, &
                tu, eigprj, pstrs, vaux, itjac1, &
                itri, iorder)
!
    if (.not. tridim) then
        eigprj(mmax, mmax) = r1
        pstrs(mmax) = stress(3)
    end if
!
! Verification of the spectral components
! ----------------------------------------
!
! Identify minor (PSTRS1) and major (PSTRS3) principal stresses
! --ok!
    ii = 1
    jj = 1
    pstrs1 = pstrs(ii)
    pstrs3 = pstrs(jj)
    do i = 2, mmax
        if (pstrs(i) .ge. pstrs1) then
            ii = i
            pstrs1 = pstrs(ii)
        end if
        if (pstrs(i) .lt. pstrs3) then
            jj = i
            pstrs3 = pstrs(jj)
        end if
    end do
    if (ii .ne. 1 .and. jj .ne. 1) mm = 1
    if (ii .ne. 2 .and. jj .ne. 2) mm = 2
    if (ii .ne. 3 .and. jj .ne. 3) mm = 3
    pstrs2 = pstrs(mm)
!
end subroutine
