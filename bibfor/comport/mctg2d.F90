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

subroutine mctg2d(stress, strain, rprops, dsidep, ii, jj, mm, &
                  edge, right, apex, outofp)
!***********************************************************************
!
!   OBJECT:
!   COMPUTE THE CONSISTENT TANGENT MATRIX FOR MOHR-COULOMB LAW
!   IN THE 2 DIMENSIONNAL CASE
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE MOHR-COULOMB
!
! IN  STRESS  : VECTEUR CONTRAINTE REACTUALISEE
! IN  STRAIN  : VECTEUR DEFORMATION A LA PREDICTION ELASTIQUE
! IN  RPROPS  : PROPIETES MECANIQUES (3)
! IN  II      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MINEURE
! IN  JJ      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MAJEURE
! IN  MM      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE INTERMEDIAIRE
! IN  EDGE    : INDICATEUR D'ACTIVATION (1) DE 2 MECANISMES
! IN  RIGHT   : INDICATEUR DE PROJECTION A DROITE (1)
!               QUAND 2 MECANISMES SONT ACITFS
! IN  APEX    : INDICATEUR D'ACTIVATION (1) DE 3 MECANISMES
! IN  OUTOFP  : COMPOSANTE HORS PLAN = | TRUE :  D_PLAN
!                                      | FALSE:  C_PLAN
!
! OUT DSIDEP  : MATRICE TANGENTE COHERENTE REACTUALISEE
!
!***********************************************************************
    implicit none
! ======================================================================
!
!
#include "asterf_types.h"
#include "asterfort/jacobi.h"
#include "asterfort/mctanp.h"
#include "asterfort/mctge2.h"
!
    real(kind=8) :: stress(6), strain(6)
    real(kind=8) :: rprops(6), dsidep(6, 6)
    real(kind=8) :: edge, right, apex
    aster_logical :: outofp
!
! Declaration of integer type variables
    integer(kind=8) :: itri, iorder, mmax, nmax, mxiter
    integer(kind=8) :: ii, jj, mm, itjac1
!
! Declaration of integer type variables
!     aster_logical :: epflag
!
    parameter(mmax=3, nmax=6)
!
! Declaration of vector and matrix type variables
    real(kind=8) :: dpstrs(mmax, mmax), pstra(mmax), pstrs(mmax)
    real(kind=8) :: eigprj(mmax, mmax), eigxpr(mmax, mmax)
    real(kind=8) :: small, vaux(mmax), tu(nmax), tr(nmax), t1(nmax)
    real(kind=8) :: r0, r1, r2, r3, r4, sqr, dmax1, tol, refe
!
! Declaration of constant variables
    data r0, r1, r2, r3, r4, small, tol, sqr/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 1.d-06, 1.d-10,&
     &    1.4142135623730951d0/
    data mxiter/50/
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
!
! Initialize unit matrix = (1 0 0 1 0 1) for Jacobi
    t1(:) = r0
    t1(1) = r1
    t1(4) = r1
    t1(6) = r1
!
    eigprj(:, :) = r0
! Matrix  TR = (SIXX SIXY SIYY) for Jacobi
! Produce EIGPRJ: Base Projection Matrix from initial base
!                 to principal directions base
!         PSTRS : principal stresses
    tr(1) = stress(1)
    tr(2) = stress(4)/sqr
    tr(3) = r0
    tr(4) = stress(2)
    tr(5) = r0
    tr(6) = r0
!
! Unit matrix = (1 0 0 1 0 1) for Jacobi
    tu(1:nmax) = t1(1:nmax)
!
    call jacobi(mmax, mxiter, tol, small, tr, &
                tu, eigprj, pstrs, vaux, itjac1, &
                itri, iorder)
!
    pstrs(3) = stress(3)
    eigprj(mmax, mmax) = r1
!
!
! Matrix  TR = (EPXX EPXY EPYY) for Jacobi
! Produce EIGXPR: Base Projection Matrix from initial base
!                 to principal directions base
!         PSTRA : principal strains
    tr(1) = strain(1)
    tr(2) = strain(4)
    tr(4) = strain(2)
    tr(3) = r0
    tr(5) = r0
    tr(6) = r0
!
! Unit matrix = (1 0 0 1 0 1) for Jacobi
    tu(1:nmax) = t1(1:nmax)
    eigxpr(:, :) = r0
!
    call jacobi(mmax, mxiter, tol, small, tr, &
                tu, eigxpr, pstra, vaux, itjac1, &
                itri, iorder)
!
    pstra(3) = strain(3)
    eigxpr(mmax, mmax) = r1
!
!
! Compute Consistent Plastic Jacobian Matrix in the eigenbasis
! Inputs: PSTRA ,RPROPS,EDGE  ,RIGHT ,APEX
! Output: DPSTRS
! ------------------------------------------------------------------
    call mctanp(dpstrs, rprops, ii, jj, mm, edge, right, apex)
!
! Check for repeated eigenvalues of strain
    refe = dmax1(abs(pstra(1)), abs(pstra(2)))*small
    if (abs(pstra(1)-pstra(2)) .lt. refe) then
        edge = r1
    else
        edge = r0
    end if
!
! Inputs:
!
!   DPSTRS : DSIGMA/DEPSI in the eigenbasis
!   PSTRA  : Strain eigenvalues vector
!   PSTRS  : Stress eigenvalues vector
!   EIGPRJ : Eigen directions Matrix (3x3)
!   EDGE   : Indicator of projection to an edge
!
! Output:
!
!   DSIDEP : Consistent Plastic Jacobian Matrix in the cartesian basis
! --------------------------------------------------------------------
    call mctge2(dpstrs, dsidep, eigxpr, pstra, pstrs, edge, outofp)
!
end subroutine
