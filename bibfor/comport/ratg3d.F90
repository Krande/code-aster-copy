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

subroutine ratg3d(stress, strain, rprops, dsidep, ii, jj, mm, &
                  edge, apex, codret)
!***********************************************************************
!
!   OBJECT:
!   COMPUTE THE CONSISTENT TANGENT MATRIX FOR MOHR-COULOMB LAW
!   IN THE 3 DIMENSIONNAL CASE
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE MOHR-COULOMB
!
! IN  STRESS  : VECTEUR CONTRAINTE REACTUALISEE
! IN  STRAIN  : VECTEUR DEFORMATION A LA PREDICTION ELASTIQUE
! IN  RPROPS  : PROPIETES MECANIQUES (3)
! IN  EDGE    : INDICATEUR D'ACTIVATION (1) DE 2 MECANISMES
! IN  APEX    : INDICATEUR D'ACTIVATION (1) DE 3 MECANISMES
! IN  II      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MINEURE
! IN  JJ      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MAJEURE
! IN  MM      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE INTERMEDIAIRE
!
! OUT DSIDEP  : MATRICE TANGENTE COHERENTE REACTUALISEE
! OUT CODRET  : CODE RETOUR
!
!***********************************************************************
    implicit none
! ======================================================================
!
    real(kind=8) :: stress(6)
    real(kind=8) :: strain(6)
    real(kind=8) :: rprops(6)
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: edge
    real(kind=8) :: apex
    integer(kind=8) :: codret
!
#include "asterf_types.h"
#include "asterfort/jacobi.h"
#include "asterfort/mcordo.h"
#include "asterfort/ratanp.h"
#include "asterfort/mctgep.h"
!
! Declaration of integer type variables
    integer(kind=8) :: itri, iorder, mmax, nmax, mxiter, itjac1
    integer(kind=8) :: ii, jj, mm
!
! Declaration of integer type variables
!     aster_logical :: epflag
!
    parameter(mmax=3, nmax=6)
!
    real(kind=8) :: dpstrs(mmax, mmax), pstra(mmax), pstrs(mmax)
    real(kind=8) :: eigprj(mmax, mmax), eigxpr(mmax, mmax)
    real(kind=8) :: small, vaux(mmax), tu(nmax), tr(nmax), t1(nmax)
    real(kind=8) :: r0, r1, r2, r3, r4, sqr, tol
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
! Matrix  TR = (SIXX SIXY SIXZ SIYY SIYZ SIZZ) for Jacobi
! Produce EIGPRJ: Base Projection Matrix from initial base
!                 to principal directions base
!         PSTRS : principal stresses
    tr(1) = stress(1)
    tr(2) = stress(4)/sqr
    tr(3) = stress(5)/sqr
    tr(4) = stress(2)
    tr(5) = stress(6)/sqr
    tr(6) = stress(3)
! Unit matrix = (1 0 0 1 0 1) for Jacobi
    tu(1:nmax) = t1(1:nmax)
!
    call jacobi(mmax, mxiter, tol, small, tr, &
                tu, eigprj, pstrs, vaux, itjac1, &
                itri, iorder)
!
! Matrix  TR = (EPXX EPXY EPXZ EPYY EPYZ EPZZ) for Jacobi
! Produce EIGXPR: Base Projection Matrix from initial base
!                 to principal directions base
!         PSTRA : principal strains
    tr(1) = strain(1)
    tr(2) = strain(4)
    tr(3) = strain(5)
    tr(4) = strain(2)
    tr(5) = strain(6)
    tr(6) = strain(3)
! Unit matrix = (1 0 0 1 0 1) for Jacobi
    tu(1:nmax) = t1(1:nmax)
!
    call jacobi(mmax, mxiter, tol, small, tr, &
                tu, eigxpr, pstra, vaux, itjac1, &
                itri, iorder)
!
! Compute Consistent Plastic Jacobian Matrix in the eigenbasis
! Inputs: PSTRA ,RPROPS,EDGE  ,RIGHT ,APEX
! Output: DPSTRS
! ------------------------------------------------------------------
    call ratanp(dpstrs, rprops, ii, jj, mm, edge, apex)
!
! Print derivative tensor
!
! Check for repeated eigenvalues of strain and re-ordering
    call mcordo(dpstrs, pstrs, pstra, eigxpr, edge, apex, codret)
!
    if (codret .eq. 1) goto 999
!
! Input Elastic Stiffness Vector in the cartesian base
! -------------------------------------------------------
    tr(1) = strain(1)
    tr(4) = strain(4)
    tr(5) = strain(5)
    tr(2) = strain(2)
    tr(6) = strain(6)
    tr(3) = strain(3)
! Inputs:
!   DPSTRS : DSIGMA/DEPSI in the eigenbasis
!   PSTRA  : Strain eigenvalues vector
!   PSTRS  : Stress eigenvalues vector
!   TR     : Total Predictor Strain vector
!   EIGPRJ : Eigen directions Matrix (3x3)
!   EDGE   : Indicator of projection to an edge
!   APEX   : Indicator of projection to the apex
!
! Output:
!   DSIDEP : Consistent Plastic Jacobian Matrix in the cartesian basis
! --------------------------------------------------------------------
    call mctgep(dpstrs, dsidep, pstra, pstrs, tr, eigxpr, edge, apex)
!
999 continue
end subroutine
